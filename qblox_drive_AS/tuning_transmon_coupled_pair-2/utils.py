"""Utility functions for executing Schedules on Qblox hardware."""

from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING

from qcodes.instrument import find_or_create_instrument

from qblox_instruments import Cluster, ClusterType
from quantify_core.measurement.control import MeasurementControl
from quantify_core.visualization.instrument_monitor import InstrumentMonitor
from quantify_core.visualization.pyqt_plotmon import PlotMonitor_pyqt as PlotMonitor
from quantify_scheduler import InstrumentCoordinator, Schedule, SerialCompiler
from quantify_scheduler.backends.types.qblox import ClusterDescription
from quantify_scheduler.instrument_coordinator.components.qblox import ClusterComponent
from quantify_scheduler.schedules.schedule import CompiledSchedule

if TYPE_CHECKING:
    from qcodes.instrument import Instrument
    from qcodes.parameters import InstrumentRefParameter
    from xarray import Dataset

    from quantify_scheduler import QuantumDevice
    from quantify_scheduler.backends.qblox_backend import QbloxHardwareCompilationConfig


DEFAULT_QUANTUM_DEVICE: QuantumDevice | None = None


def initialize_hardware(
    quantum_device: QuantumDevice | None = None, ip: str | None = None, live_plotting: bool = False
) -> tuple[MeasurementControl, InstrumentCoordinator, Cluster]:
    """
    Initialize MeasurementControl and InstrumentCoordinator from QuantumDevice.

    Parameters
    ----------
    quantum_device : QuantumDevice | None, optional
        target quantum device, by default None
    ip : str | None, optional
        ip address of the qblox cluster. Will use a dummy cluster if None, by default None
    live_plotting : bool, optional
        wether live plotting should be enabled, by default False

    Returns
    -------
    tuple[MeasurementControl, InstrumentCoordinator]

    Raises
    ------
    ValueError
        Neither QuantumDevice nor global default are provided.

    """
    if quantum_device is None:
        if DEFAULT_QUANTUM_DEVICE is None:
            raise ValueError("Either provide a QuantumDevice or set the global default")
        else:
            quantum_device = DEFAULT_QUANTUM_DEVICE
    config: QbloxHardwareCompilationConfig = quantum_device.generate_hardware_compilation_config()
    all_instruments = []
    for name, instrument in config.hardware_description.items():
        if isinstance(instrument, ClusterDescription):
            dummy_config = {}
            dummy_type_by_str = {
                "QCM": ClusterType.CLUSTER_QCM,
                "QRM": ClusterType.CLUSTER_QRM,
                "QCM_RF": ClusterType.CLUSTER_QCM_RF,
                "QRM_RF": ClusterType.CLUSTER_QRM_RF,
                "QTM": ClusterType.CLUSTER_QTM,
                "QDM": ClusterType.CLUSTER_QDM,
            }
            instrument_ip = ip if ip is not None else getattr(instrument, "ip", None)
            if instrument_ip is None:
                for slot, module in instrument.modules.items():
                    dummy_config[int(slot)] = dummy_type_by_str[module.instrument_type]

            cluster = find_or_create_instrument(
                Cluster,
                recreate=True,
                name=name,
                identifier=instrument_ip,
                dummy_cfg=dummy_config,
            )
            all_instruments.append(cluster)

    # Associate measurement controller and instrument coordinator with the quantum device
    def instantiate_instrument(
        parameter: InstrumentRefParameter, instr_type: Instrument, name: str
    ) -> Instrument:
        try:
            inst = parameter.get_instr()
        except KeyError as e:
            if "Instrument with name " in str(e):
                inst = find_or_create_instrument(instr_type, recreate=True, name=name)
            else:
                raise e
            parameter(name)
        return inst

    meas_ctrl = instantiate_instrument(
        quantum_device.instr_measurement_control, MeasurementControl, "meas_ctrl"
    )
    ic = instantiate_instrument(
        quantum_device.instr_instrument_coordinator, InstrumentCoordinator, "ic"
    )

    for cluster in all_instruments:
        # Add cluster to instrument coordinator
        ic_cluster = ClusterComponent(cluster)
        with contextlib.suppress(ValueError):
            ic.add_component(ic_cluster)

    if live_plotting:
        instantiate_instrument(meas_ctrl.instr_plotmon, PlotMonitor, "PlotMonitor")
        find_or_create_instrument(InstrumentMonitor, recreate=True, name="Instrument_Monitor")

    return (meas_ctrl, ic, cluster)


def run(schedule: Schedule, quantum_device: QuantumDevice) -> Dataset:
    if not isinstance(schedule, CompiledSchedule):
        compiler = SerialCompiler(name="auto_compiler", quantum_device=quantum_device)
        schedule = compiler.compile(schedule=schedule)
    ic: InstrumentCoordinator = quantum_device.instr_instrument_coordinator.get_instr()
    ic.prepare(schedule)
    ic.start()
    return ic.retrieve_acquisition()
