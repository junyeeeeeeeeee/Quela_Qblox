from __future__ import annotations
from typing import Optional, Union, List
import warnings

import numpy as np

from numpy.typing import NDArray

from quantify_scheduler.operations.operation import Operation
from quantify_scheduler.resources import BasebandClockResource
from quantify_scheduler.operations.pulse_library import ReferenceMagnitude

def Sin(
        t: Union[np.ndarray, List[float]],
        amp: Union[float, complex],
) -> np.ndarray:
    
    # waveform = amp * 1.0j * np.exp(1.0j * t/t[-1] * np.pi)
    waveform = -amp * np.sin(t/t[-1] * np.pi)
    return waveform

class SinPulse(Operation):  # pylint: disable=too-many-ancestors
    """
    Sin

    Parameters
    ----------
    amp
        amp
    duration
        Duration of the pulse in seconds.
    port
        Port of the pulse.
    clock
        Clock used to modulate the pulse.
    reference_magnitude
        Scaling value and unit for the unitless amplitude. Uses settings in
        hardware config if not provided.
    t0
        Time in seconds when to start the pulses relative to the start time
        of the Operation in the Schedule.
    """

    def __init__(
        self,
        amp: float,
        duration: float,
        phase: float,
        port: str,
        clock: str = BasebandClockResource.IDENTITY,
        reference_magnitude: Optional[ReferenceMagnitude] = None,
        t0: float = 0,
    ):
        super().__init__(name=self.__class__.__name__)
        self.data["pulse_info"] = [
            {
                "wf_func": "Modularize.support.custumized_pulse_library.Sin",
                "amp": amp,
                "reference_magnitude": reference_magnitude,
                "duration": duration,
                "phase": phase,
                "t0": t0,
                "clock": clock,
                "port": port,
            }
        ]
        self._update()

    def __str__(self) -> str:
        pulse_info = self.data["pulse_info"][0]
        return self._get_signature(pulse_info)