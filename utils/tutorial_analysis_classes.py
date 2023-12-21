from typing import Union

import lmfit
import matplotlib.pyplot as plt
import numpy as np
import quantify_core.data.handling as dh
from numpy.typing import NDArray
from quantify_core.analysis import base_analysis as ba
from quantify_core.analysis import fitting_models as fm
from quantify_core.visualization import mpl_plotting as qpl
from quantify_core.visualization.SI_utilities import format_value_string
from tqdm.auto import tqdm


def _helper_sinus(
    x: Union[NDArray, float], *, w: float, amp: float, phs: float, offset: float
) -> Union[NDArray, float]:
    """Return values of a parameterized sinusoid."""
    return float(amp) * np.sin(x * float(w) + float(phs)) + float(offset)


def _helper_find_zero(*, n: int, phs: float, w: float) -> float:
    """Return the n'th zero of the derivative of `_helper_sinus()`."""
    return (-2.0 * phs - 2 * np.pi * int(n) + np.pi) / (2.0 * w)


class ResonatorFluxSpectroscopyAnalysis(ba.BaseAnalysis):
    """
    Analysis class for resonator flux spectroscopy.

    ((y0, y1),(x0, x1)) maps to ((mag, phase), (freqs, offsets)).
    """

    def _find_resonance_points_naive(self) -> NDArray:
        """Extract resonance frequencies with arg minima of magnitude image."""
        # Fetch setpoint data
        frequencies = self.dataset_processed.x0

        # Fetch magnitude image data
        magnitude = self.dataset_processed.y0

        # Guess that the resonance is the lowest point in magnitude
        # (leftmost when multiple)
        argmin_s21 = np.argmin(magnitude.data, axis=0)
        return frequencies.data[argmin_s21]

    def _find_resonance_points_precise(self) -> NDArray:
        """Extract resonance frequencies with resonator models."""
        # Convert magnitude, phase data to complex point data
        S21 = self.dataset_processed.y0.data * np.cos(
            np.deg2rad(self.dataset_processed.y1.data)
        ) + 1j * self.dataset_processed.y0.data * np.sin(
            np.deg2rad(self.dataset_processed.y1.data)
        )
        S21 = np.transpose(S21)

        # Fetch setpoint data
        frequencies = self.dataset_processed.x0

        # Instantiate resonator models
        resonators = [fm.ResonatorModel() for _ in range(S21.shape[0])]
        fitted_params = []
        fit_results_all = []

        for offs_idx in tqdm(range(S21.shape[0]), desc="Fitting resonator models"):
            spectroscopy = S21[offs_idx, :]
            res = resonators[offs_idx]

            # Use last found parameters as best guess
            # (this is significantly faster than re-guessing every time)
            if len(fitted_params) > 0:
                guess = fitted_params[-1]

            # First time guess, use the guess from ResonatorModel
            else:
                guess = res.guess(spectroscopy, f=frequencies.data)

            fit_results = res.fit(spectroscopy, f=frequencies.data, params=guess)
            fitted_params.append(fit_results.params)
            fit_results_all.append(fit_results)

        return (np.asarray([fp["fr"].value for fp in fitted_params]), fit_results_all)

    def run(self, fit_method: str = "fast", sweetspot_index: int = 0) -> object:
        """Execute all analysis steps of this class and return it."""
        # Execute the analysis steps directly
        self.process_data()
        self.run_fitting(fit_method=fit_method)
        self.analyze_fit_results(sweetspot_index=sweetspot_index)
        self.create_figures()
        self.adjust_figures()
        self.save_figures()
        self.save_quantities_of_interest()
        self.save_processed_dataset()
        self.save_fit_results()

        # Return the analysis object
        return self

    def process_data(self) -> None:
        """Process the data so that the analysis can make assumptions on the format."""
        if self.dataset.attrs["grid_2d"]:
            self.dataset_processed = dh.to_gridded_dataset(self.dataset)

    def run_fitting(self, *, verbose: bool = False, fit_method: str = "fast") -> None:
        """Fits a sinusoidal model to the frequency response vs. flux offset."""
        # Fetch setpoint data
        offsets = self.dataset_processed.x1

        # FIXME: Workaround
        self._fit_method = str.lower(fit_method)

        # Find resonance frequencies for each spectroscopy
        if str.lower(fit_method) == "precise":
            f_minima, fit_results_all = self._find_resonance_points_precise()

            self.fit_results.update(
                {f"resonator_model_{r}": frr for r, frr in enumerate(fit_results_all)}
            )

        else:
            f_minima = self._find_resonance_points_naive()

        # Center the resonance point data for model fitting
        f_minima_mean = f_minima.mean()
        f_minima_centered = f_minima - f_minima_mean

        # Fit a sinusoidal model to centered data
        sine_model = lmfit.models.SineModel()
        guess = sine_model.guess(f_minima_centered, x=offsets.data)
        result = sine_model.fit(f_minima_centered, x=offsets.data, params=guess)

        # Also add the offset with estimate standard error
        result.params.add("offset", value=f_minima_mean, vary=False)
        result.params["offset"].stderr = f_minima.std() / np.sqrt(f_minima.shape[0])

        # Verbosity:
        if verbose:
            print(result.fit_report() + "\n")

        self.fit_results.update({"sin": result})

    def analyze_fit_results(self, sweetspot_index: int) -> None:
        """Check the fit success and populate :code:`.quantities_of_interest`."""
        self.quantities_of_interest = {}

        # Set fitted sinus parameters as quantities of interest
        for parameter, value in self.fit_results["sin"].params.items():
            self.quantities_of_interest[parameter] = ba.lmfit_par_to_ufloat(value)

        # Find the first zero of the derivative of the fitted sinusoidal model
        self.quantities_of_interest["offset_0"] = _helper_find_zero(
            n=sweetspot_index,
            w=max(1e-30, abs(self.quantities_of_interest["frequency"])),
            phs=self.quantities_of_interest["shift"],
        )

        # If you have fitted resonators already, take closest to the optimal flux.
        if self._fit_method == "precise":
            closest_resonator_idx = np.abs(
                self.dataset_processed.x1.data
                - self.quantities_of_interest["offset_0"].nominal_value
            ).argmin()
            self.quantities_of_interest[
                "sweetspot_resonator_idx"
            ] = closest_resonator_idx
            self.fit_results["sweetspot_resonator_model"] = self.fit_results[
                f"resonator_model_{closest_resonator_idx}"
            ]

        # Evaluate the model at that point to extract the corresponding frequency
        self.quantities_of_interest["freq_0"] = float(
            _helper_sinus(
                x=self.quantities_of_interest["offset_0"].nominal_value,  # type: ignore
                w=self.quantities_of_interest["frequency"].nominal_value,
                amp=self.quantities_of_interest["amplitude"].nominal_value,
                phs=self.quantities_of_interest["shift"].nominal_value,
                offset=self.quantities_of_interest["offset"].nominal_value,
            )
        )

        # If there is a problem with the fit, display an error message in the text box.
        # Otherwise, display the parameters as normal.
        fit_warning = ba.wrap_text(ba.check_lmfit(self.fit_results["sin"]))
        if fit_warning is None:
            self.quantities_of_interest["fit_success"] = True

            text_msg = "Summary:\n"
            text_msg += format_value_string(
                "f_0",
                self.quantities_of_interest["freq_0"],
                unit="Hz",
                end_char="\n",
            )
            text_msg += format_value_string(
                "offset_0",
                self.quantities_of_interest["offset_0"],
                unit=self.dataset_processed.x1.units,
                end_char="\n",
            )

        else:
            text_msg = ba.wrap_text(fit_warning)
            self.quantities_of_interest["fit_success"] = False

        self.quantities_of_interest["fit_msg"] = text_msg

    def create_figures(self) -> None:
        """Generate plot of magnitude and phase images, with superposed model fit."""
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=plt.figaspect(1 / 2))

        # Plot using xarrays plotting function
        self.dataset_processed.y0.plot(ax=ax0)  # type: ignore
        np.mod(self.dataset_processed.y1, 360.0).plot(ax=ax1)  # type: ignore

        # Interpolate with the model for nice curve
        interp_offsets = np.arange(
            start=self.dataset_processed.x1.data.min(),
            stop=self.dataset_processed.x1.data.max(),
            step=np.diff(self.dataset_processed.x1)[0] / 2.0,
        )

        for ax in (ax0, ax1):
            # Plot model interpolation
            ax.plot(
                interp_offsets,
                _helper_sinus(
                    x=interp_offsets,
                    w=self.quantities_of_interest["frequency"].nominal_value,
                    amp=self.quantities_of_interest["amplitude"].nominal_value,
                    phs=self.quantities_of_interest["shift"].nominal_value,
                    offset=self.quantities_of_interest["offset"].nominal_value,
                ),
                color="red",
            )

            # Plot also the zero which we found
            ax.axvline(
                self.quantities_of_interest["offset_0"].nominal_value,
                color="red",
                ls="--",
            )

        # Put a text box summarizing the QOI in the middle of the figure
        qpl.plot_textbox(
            ax=ax1,
            text=self.quantities_of_interest["fit_msg"],
            x=1.34,  # make some space for colorbar
        )

        qpl.set_suptitle_from_dataset(fig, self.dataset)  # type: ignore
        ax0.set_title("Magnitude")
        ax1.set_title("Phase")
        fig.tight_layout()  # type: ignore

        self.figs_mpl["rfs"] = fig  # type: ignore
        self.axs_mpl["rfs"] = (  # type: ignore
            ax0,
            ax1,
        )


class QubitFluxSpectroscopyAnalysis(ba.BaseAnalysis):
    """
    Analysis class for qubit flux spectroscopy.

    ((y0, y1),(x0, x1)) maps to ((mag, phase), (freqs, offsets)).
    """

    def _to_transmission(self) -> NDArray:
        """Return S_21 of magnitude and phase data."""
        S21 = self.dataset_processed.y0.data * np.cos(
            np.deg2rad(self.dataset_processed.y1.data)
        ) + 1j * self.dataset_processed.y0.data * np.sin(
            np.deg2rad(self.dataset_processed.y1.data)
        )
        return np.transpose(S21)

    def _find_frequencies(self, offsets: NDArray, frequencies: NDArray) -> NDArray:
        """Extract offsets and frequencies for fitting the quadratic model."""
        s21 = self._to_transmission()

        # output storage
        pts = np.empty((offsets.shape[0], 2))
        pts[:] = np.nan

        for offset_idx, row in enumerate(s21):
            row_mean = row.mean()

            # outlier if larger than tau away from the mean
            row_tau = 3 * row.std()
            delta = np.sqrt(
                (row_mean.real - row.real) ** 2 + (row_mean.imag - row.imag) ** 2
            )

            # take frequencies of outliers for every offset
            f = frequencies[delta > row_tau]

            if len(f) > 0:
                # copy independent variable (offset)
                pts[offset_idx, 0] = offsets[offset_idx]

                # dependent variable
                # need to fit as a function, so merge datapoints with expected value
                pts[offset_idx, 1] = f.mean()

            # if no outliers are found, row in output will be [np.nan,np.nan]

        # filter to remove data points where there is no peak data
        return pts[~np.isnan(pts).any(axis=1)]

    def run(self, /, *, verbose: bool = False) -> object:
        """Execute all analysis steps of this class and return it."""
        # Execute the analysis steps directly
        self.process_data()
        self.run_fitting(verbose=verbose)
        self.analyze_fit_results()
        self.create_figures()
        self.adjust_figures()
        self.save_figures()
        self.save_quantities_of_interest()
        self.save_processed_dataset()
        self.save_fit_results()

        # Return the analysis object
        return self

    def process_data(self) -> None:
        """Process the data so that the analysis can make assumptions on the format."""
        if self.dataset.attrs["grid_2d"]:
            self.dataset_processed = dh.to_gridded_dataset(self.dataset)

    def run_fitting(self, *, verbose: bool = False) -> None:
        """Fits a quadratic model to the frequency response vs. flux offset."""
        u = self._find_frequencies(
            frequencies=self.dataset_processed.x0.data,
            offsets=self.dataset_processed.x1.data,
        )

        # Fit a model to data
        quad_model = lmfit.models.QuadraticModel()
        guess = quad_model.guess(u[:, 1], x=u[:, 0])
        result = quad_model.fit(u[:, 1], x=u[:, 0], params=guess)

        # Verbosity:
        if verbose:
            print(result.fit_report() + "\n")

        self.fit_results.update({"poly2": result})

    def analyze_fit_results(self) -> None:
        """Check the fit success and populate :code:`.quantities_of_interest`."""
        self.quantities_of_interest = {}

        # Set fitted parameters as quantities of interest
        for parameter, value in self.fit_results["poly2"].params.items():
            self.quantities_of_interest[parameter] = ba.lmfit_par_to_ufloat(value)

        # If there is a problem with the fit, display an error message in the text box.
        # Otherwise, display the parameters as normal.
        fit_warning = ba.wrap_text(ba.check_lmfit(self.fit_results["poly2"]))
        if fit_warning is None:
            self.quantities_of_interest["fit_success"] = True

            text_msg = "Summary:\n"

            a = self.quantities_of_interest["a"]
            b = self.quantities_of_interest["b"]
            c = self.quantities_of_interest["c"]

            off_0_unc = -b / (2.0 * a)
            frq_0_unc = a * (off_0_unc**2) + b * off_0_unc + c

            self.quantities_of_interest["offset_0"] = off_0_unc
            self.quantities_of_interest["freq_0"] = frq_0_unc

            text_msg += format_value_string(
                "offset_0",
                off_0_unc,
                unit=self.dataset_processed.x1.units,
                end_char="\n",
            )

            text_msg += format_value_string(
                "freq_0",
                frq_0_unc,
                unit="Hz",
                end_char="\n",
            )

        else:
            text_msg = ba.wrap_text(fit_warning)
            self.quantities_of_interest["fit_success"] = False

        self.quantities_of_interest["fit_msg"] = text_msg

    def create_figures(self) -> None:
        """Generate plot of magnitude and phase images, with superposed model fit."""
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=plt.figaspect(1 / 2.8))

        # Plot using xarrays plotting function
        self.dataset_processed.y0.plot(ax=ax0)  # type: ignore
        np.mod(self.dataset_processed.y1, 360.0).plot(ax=ax1)  # type: ignore

        a = self.quantities_of_interest["a"].nominal_value
        b = self.quantities_of_interest["b"].nominal_value
        c = self.quantities_of_interest["c"].nominal_value

        interp_offsets = np.arange(
            start=self.dataset_processed.x1.data.min(),
            stop=self.dataset_processed.x1.data.max(),
            step=np.diff(self.dataset_processed.x1.data)[0] / 2.0,
        )
        for ax in (ax0, ax1):
            ax.plot(
                interp_offsets,
                a * (interp_offsets**2) + b * interp_offsets + c,
                color="red",
                ls="-",
            )
            ax.axvline(
                self.quantities_of_interest["offset_0"].nominal_value,
                color="red",
                ls="--",
            )

        # Put a text box summarizing the QOI in the middle of the figure
        qpl.plot_textbox(
            ax=ax1,
            text=self.quantities_of_interest["fit_msg"],
            x=1.34,  # make some space for colorbar
        )

        qpl.set_suptitle_from_dataset(fig, self.dataset)  # type: ignore
        ax0.set_title("Magnitude")
        ax1.set_title("Phase")
        fig.tight_layout()  # type: ignore

        self.figs_mpl["qfs"] = fig  # type: ignore
        self.axs_mpl["qfs"] = (  # type: ignore
            ax0,
            ax1,
        )


class ConditionalOscillationAnalysis(ba.BaseAnalysis):
    """
    Perform the analysis of a conditional oscillation experiment.

    For a reference to the conditional oscillation experiment, please
    see section D in the supplemental material of
    this paper: https://arxiv.org/abs/1903.02492
    """

    # pylint: disable=no-member
    # pylint: disable=invalid-name

    def process_data(self) -> None:
        """Process the data so that the analysis can make assumptions on the format."""
        # self.dataset_processed = dh.to_gridded_dataset(self.dataset)
        self.dataset_processed = self.dataset

        # magnitude
        assert self.dataset_processed.y0.units == "V"  # off
        assert self.dataset_processed.y2.units == "V"  # off
        assert self.dataset_processed.y4.units == "V"  # on
        assert self.dataset_processed.y6.units == "V"  # on

        # phase
        assert self.dataset_processed.y1.units == "deg"  # off
        assert self.dataset_processed.y3.units == "deg"  # off
        assert self.dataset_processed.y5.units == "deg"  # on
        assert self.dataset_processed.y7.units == "deg"  # on

    def run(self) -> object:
        """Execute all analysis steps of this class and return it."""
        # Execute the analysis steps directly
        self.process_data()
        self.run_fitting()
        self.analyze_fit_results()
        self.create_figures()
        self.adjust_figures()
        self.save_figures()
        self.save_quantities_of_interest()
        self.save_processed_dataset()
        self.save_fit_results()

        # Return the analysis object
        return self

    def run_fitting(self) -> None:
        """Fit model to the data."""
        self.fit_results = {}

        def _center_and_fit_sinus(y: NDArray, x: NDArray) -> lmfit.model.ModelResult:
            # Center the data
            y_centered = y - y.mean()

            # Fit a sinusoidal model to centered data (ON)
            sine_model = lmfit.models.SineModel()
            guess = sine_model.guess(y_centered, x=x)
            result = sine_model.fit(y_centered, x=x, params=guess)

            return result

        self.fit_results.update(
            {
                "sin_off": _center_and_fit_sinus(
                    self.dataset_processed.y0.data, self.dataset_processed.x0.data
                )
            }
        )

        self.fit_results.update(
            {
                "sin_on": _center_and_fit_sinus(
                    self.dataset_processed.y4.data, self.dataset_processed.x0.data
                )
            }
        )

    def analyze_fit_results(self) -> None:
        """
        Check fit success and populates :code:`.quantities_of_interest`.

        "leak" estimate measure of the leakage to the |2> state
        for the HF qubit.

        "phi_2q_deg" conditional phase of the 2-qubit gate.

        'off_amplitude' amplitude of the sinusoidal LF qubit measurement
        in the ON variant.

        'off_frequency' frequency of the sinusoidal LF qubit measurement
        in the ON variant.

        'off_shift' phase of the sinusoidal LF qubit measurement in
        the ON variant.

        'off_hf_level' average value of the HF qubit measurement in the
        ON variant.

        'off_offset' constant offset of the sinusoidal LF qubit measurement
        in the ON variant.

        'on_amplitude' amplitude of the sinusoidal LF qubit measurement
        in the OFF variant.

        'on_frequency' frequency of the sinusoidal LF qubit measurement
        in the OFF variant.

        'on_shift'  phase of the sinusoidal LF qubit measurement in
        the OFF variant.

        'on_hf_level' average value of the HF qubit measurement in the
        OFF variant.

        'on_offset' constant offset of the sinusoidal LF qubit measurement
        in the OFF variant.

        """
        self.quantities_of_interest = {}

        def _add_center(
            param_name: str, data: NDArray, params: lmfit.parameter.Parameters
        ) -> None:
            params.add(param_name, value=data.mean(), vary=False)
            params[param_name].stderr = data.std() / np.sqrt(data.shape[0])

        def _store_params(data: NDArray, mode: str) -> None:
            result = self.fit_results["sin_" + str(mode)]

            _add_center("offset", data, result.params)

            # Store all parameters as quantities of interest
            for p, value in result.params.items():
                self.quantities_of_interest[
                    str(mode) + "_" + str(p)
                ] = ba.lmfit_par_to_ufloat(value)

        # Extract proper leakage estimator (TODO: Relate quantity to |1> proportion)
        # store in parameters struct of model
        _add_center(
            "hf_level",
            self.dataset_processed.y6.data,
            self.fit_results["sin_on"].params,
        )
        _add_center(
            "hf_level",
            self.dataset_processed.y2.data,
            self.fit_results["sin_off"].params,
        )

        _store_params(self.dataset_processed.y0.data, mode="off")
        _store_params(self.dataset_processed.y4.data, mode="on")

        # Extract conditional phase in degrees
        self.quantities_of_interest["phi_2q_deg"] = (
            self.quantities_of_interest["on_shift"]
            - self.quantities_of_interest["off_shift"]
        )
        self.quantities_of_interest["phi_2q_deg"] *= 180.0 / np.pi  # convert to degrees
        self.quantities_of_interest["phi_2q_deg"] = np.mod(
            self.quantities_of_interest["phi_2q_deg"], 360.0
        )

        self.quantities_of_interest["leak"] = ba.lmfit_par_to_ufloat(
            self.fit_results["sin_on"].params["hf_level"]
        ) - ba.lmfit_par_to_ufloat(
            self.fit_results["sin_off"].params["hf_level"]
        )  # type: ignore

        # Make fit message
        text_msg = "Summary:\n\n"
        text_msg += format_value_string(
            "Phi",
            self.quantities_of_interest["phi_2q_deg"],
            unit="deg",
            end_char="\n",
        )
        text_msg += format_value_string(
            "Leakage",
            self.quantities_of_interest["leak"],
            unit="V",
            end_char="\n",
        )
        self.quantities_of_interest["fit_msg"] = text_msg

    def create_figures(self) -> None:
        """
        Generate figures of interest.

        matplolib figures and axes objects are added to
        the .figs_mpl and .axs_mpl dictionaries., respectively.
        """
        fig, axs = plt.subplots(1, 2, figsize=plt.figaspect(1 / 2), sharey=True)

        # plot lf measurement for conditional phase
        self.dataset_processed.y0.plot(
            ax=axs[0], marker=".", label="low freq. qb off", color="C0"
        )
        self.dataset_processed.y4.plot(
            ax=axs[0], marker=".", label="low freq. qb on", color="C1"
        )

        # plot hf measurement for leakage estimator
        self.dataset_processed.y2.plot(
            ax=axs[1], marker=".", label="high freq. qb off", color="C0"
        )
        self.dataset_processed.y6.plot(
            ax=axs[1], marker=".", label="high freq. qb on", color="C1"
        )

        # Interpolate with the model for nice curve
        interp_x = np.arange(
            start=self.dataset_processed.x0.data.min(),
            stop=self.dataset_processed.x0.data.max(),
            step=np.diff(self.dataset_processed.x0)[0] / 2.0,
        )

        for var in ("on", "off"):
            axs[0].plot(
                interp_x,
                _helper_sinus(
                    x=interp_x,
                    w=self.quantities_of_interest[var + "_frequency"].nominal_value,
                    amp=self.quantities_of_interest[var + "_amplitude"].nominal_value,
                    phs=self.quantities_of_interest[var + "_shift"].nominal_value,
                    offset=self.quantities_of_interest[var + "_offset"].nominal_value,
                ),
                color="red",
                ls="-",
            )

        axs[1].axhline(self.dataset_processed.y2.data.mean(), color="red", ls="--")
        axs[1].axhline(self.dataset_processed.y6.data.mean(), color="red", ls="--")

        for ax in axs:
            ax.grid(alpha=1 / 25, color="black")
            ax.legend(loc="upper left")

        qpl.set_suptitle_from_dataset(fig, self.dataset)  # type: ignore
        fig.tight_layout()  # type: ignore

        # Put a text box summarizing the QOI in the middle of the figure
        qpl.plot_textbox(
            ax=axs[1],
            text=self.quantities_of_interest["fit_msg"],
            # x=1.34, # make some space for colorbar
        )

        self.figs_mpl["ConditionalOscillationAnalysis"] = fig  # type: ignore
        self.axs_mpl["ConditionalOscillationAnalysis"] = axs  # type: ignore
