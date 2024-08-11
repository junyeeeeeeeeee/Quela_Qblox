from typing import Optional

import pandas as pd
from quantify_scheduler.device_under_test.transmon_element import BasicTransmonElement
from quantify_scheduler.helpers.collections import find_port_clock_path


def show_args(my_argument_dict: dict, title: Optional[str] = None) -> None:
    max_len_keys = max(set(map(len, my_argument_dict.keys())))
    max_len_values = max(set(map(lambda v: len(str(v)), my_argument_dict.values())))

    if title is not None:
        print(str(title))
        print("=" * (max_len_keys + max_len_values))

    print(
        "\n".join(
            [
                f"{k}\t= {v}".expandtabs(max_len_keys + 1)
                for k, v in my_argument_dict.items()
            ]
        )
    )


def show_parameters(*basic_transmon_elements: list) -> None:
    pd.set_option("display.precision", 10)

    df = pd.DataFrame()

    def _units(container: dict) -> list:
        return [f"({v.unit})" if v.unit != "" else "" for v in container.values()]

    def _values(container: dict) -> list:
        return list(map(lambda v: v(), container.values()))

    # No elements provided, just print empty table
    if not len(basic_transmon_elements):
        print(pd.DataFrame())
        return

    # Test that all provided elements are correct type
    if not all(
        [isinstance(bte, BasicTransmonElement) for bte in basic_transmon_elements]
    ):
        raise TypeError("All provided objects should be BasicTransmonElements.")

    # At least one BasicTransmonElement was provided, use this to get keys
    r = basic_transmon_elements[0]
    measure_keys = list(r.measure.parameters.keys())
    rxy_keys = list(r.rxy.parameters.keys())
    reset_keys = list(r.reset.parameters.keys())
    clock_freqs_keys = list(r.clock_freqs.parameters.keys())
    ports_keys = list(r.ports.parameters.keys())

    df["Parameter"] = (
        measure_keys + reset_keys + clock_freqs_keys + rxy_keys + ports_keys
    )
    df["Type"] = (
        ["measure"] * len(measure_keys)
        + ["reset"] * len(reset_keys)
        + ["clock_freqs"] * len(clock_freqs_keys)
        + ["rxy"] * len(rxy_keys)
        + ["ports"] * len(ports_keys)
    )
    df["Unit"] = (
        _units(r.measure.parameters)
        + _units(r.reset.parameters)
        + _units(r.clock_freqs.parameters)
        + _units(r.rxy.parameters)
        + [""] * len(ports_keys)
    )

    for bte in basic_transmon_elements:
        df[bte.name] = (
            _values(bte.measure.parameters)
            + _values(bte.reset.parameters)
            + _values(bte.clock_freqs.parameters)
            + _values(bte.rxy.parameters)
            + _values(bte.ports.parameters)
        )

    df.set_index("Parameter", inplace=True)
    print(df)


def show_readout_args(qubit, /) -> None:
    show_args(
        {k: v() for k, v in qubit.measure.parameters.items()},
        title=f"{qubit.name}.measure",
    )
    show_args(
        {k: v() for k, v in qubit.reset.parameters.items()},
        title=f"\n{qubit.name}.reset",
    )
    show_args(
        {"readout": qubit.clock_freqs.readout()}, title=f"\n{qubit.name}.clock_freqs"
    )


def show_drive_args(qubit, /) -> None:
    show_args(
        {k: v() for k, v in qubit.rxy.parameters.items()}, title=f"{qubit.name}.rxy"
    )
    show_args({"f01": qubit.clock_freqs.f01()}, title=f"\n{qubit.name}.clock_freqs")


def set_readout_attenuation(quantum_device, qubit, /, *, out_att: int, in_att: int):
    """
    Set output and input attenuation of QRM RF.

    TODO: Double check that in_att > 0 dB works on hardware.
    """
    hw_config = quantum_device.hardware_config()

    output_path = find_port_clock_path(
        hw_config, port='q:res', clock=qubit.name + ".ro"  #qubit.ports.readout()
    )

    cluster_key, module_key, output_key, _, _ = tuple(output_path)
    readout_module = hw_config[cluster_key][module_key]
    readout_module[output_key]["output_att"] = out_att
    readout_module[output_key]["input_att"] = in_att

    quantum_device.hardware_config(hw_config)


def set_drive_attenuation(quantum_device, qubit, /, *, out_att: int):
    """
    Set output attenuation of QCM RF.
    """
    hw_config = quantum_device.hardware_config()

    output_path = find_port_clock_path(
        hw_config, port=qubit.ports.microwave(), clock=qubit.name + ".01"
    )

    cluster_key, module_key, output_key, _, _ = tuple(output_path)
    drive_module = hw_config[cluster_key][module_key]
    drive_module[output_key]["output_att"] = out_att

    quantum_device.hardware_config(hw_config)
