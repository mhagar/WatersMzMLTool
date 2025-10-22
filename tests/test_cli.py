from dataclasses import dataclass

@dataclass
class MockArguments:
    input: str
    experiment_type: str
    add_ms_levels: bool
    remove_calibration_scans: bool
    centroid: bool
    calibrants: str = None
    formulae: bool = False
    calibrant_window: float = 0.1
    calibrant_min_intensity: float = 1e3
    show_calibration_plot: bool = False
    output: str = None
    gui: bool = False


def test_process_mzml_file_cli(test_file_path):
    from waters_mzml_tools import process_mzml_file_cli

    args = MockArguments(
        input=str(test_file_path),
        experiment_type='dia',
        add_ms_levels=True,
        remove_calibration_scans=True,
        centroid=True,
        calibrants="/home/mh/PycharmProjects/WatersMzMLTool/calibration_series_leuenk_glufib_fragment.txt",
    )

    process_mzml_file_cli(
        input_path=test_file_path,
        args=args  # type: ignore
    )