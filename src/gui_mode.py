"""
This script manages the GUI mode, which is launched if no
arguments are provided to the CLI mode, or if --gui is passed
"""
import pyqtgraph as pg
from PyQt5 import QtWidgets
from PyQt5.QtCore import QThread

import src.calibration as calibration
from src.gui.main_window_layout import Ui_MainWindow
from src.qt_worker import ProcessingWorker

from typing import Literal, Optional
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Parameters:
    """
    Object formed from the GUI when the user hits 'Execute'
    """
    input_files: list[Path]
    experiment_type: Literal['dia', 'dda']
    output_dir: Path
    remove_lockmass_scans: bool
    fix_ms_levels: bool
    calibrate: bool
    calibrant_file_dir: Path
    calibrant_file_formulae: bool
    calibrant_search_window: float
    calibrant_min_intsy: float


    def format_as_cli_arguments(self) -> list[str]:
        args: list[str] = []
        args += ['--experiment-type', str(self.experiment_type)]

        if self.fix_ms_levels: args += ['--add-ms-levels']
        if self.output_dir: args += ['--output', str(self.output_dir)]

        if self.calibrate:
            args += ['--calibrants', str(self.calibrant_file_dir)]
            if self.calibrant_file_formulae: args += ['--formulae']
            args += ['--calibrant-window', str(self.calibrant_search_window)]
            args += ['--calibrant-min-intensity', str(self.calibrant_min_intsy)]

        if self.remove_lockmass_scans:
            args += ['--remove-calibration-scans']

        return args


class MainWindow(
    QtWidgets.QMainWindow,
    Ui_MainWindow,
):
    def __init__(
        self
    ):
        super().__init__()
        self.setupUi(MainWindow=self)
        self.parameters: Optional[Parameters] = None
        self.worker: Optional[ProcessingWorker] = None
        self.thread: Optional[QThread] = None
        self.progress_dialog: Optional[QtWidgets.QProgressDialog] = None

    def btnAddFilePressed(self):
        filepaths, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self,
            caption='Select .mzML files for manipulation',
            filter='LC/MS data (*.mzML)',
        )

        if not filepaths:
            # User cancelled
            return

        currently_displayed_paths: list[str] = self.get_input_files()

        for filepath in filepaths:
            if not Path(filepath).exists():
                print(
                    f"{filepath} does not exist!"
                )

            if filepath not in currently_displayed_paths:
                self.listInputFiles.addItem(
                    filepath
                )

    def btnRemoveFilePressed(self):
        """
        If there's an item selected in the list widget, removes it.
        :return:
        """
        for item in self.listInputFiles.selectedItems():
            self.listInputFiles.takeItem(
                self.listInputFiles.row(item)
            )

    def btnOutputDirPressed(self):
        output_dir = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            caption="Select directory to save output .mzML files"
        )

        if output_dir:
            self.lineOutputDir.setText(
                output_dir
            )

    def btnCalibrantDirPressed(self):
        calibrant_path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            caption="Select file listing calibrant masses",
            filter="Calibrant file (*.txt)",
        )

        if calibrant_path:
            self.lineCalibrantDir.setText(
                calibrant_path
            )

    def retrieve_parameters_from_ui(self):
        self.parameters = Parameters(
            input_files=self.get_input_files(),
            experiment_type=self.comboExperimentType.currentText().lower(),
            output_dir=self.lineOutputDir.text(),
            remove_lockmass_scans=self.checkRemoveLockmassScans.isChecked(),
            fix_ms_levels=self.checkFixMSLevels.isChecked(),
            calibrate=self.groupCalibrate.isChecked(),
            calibrant_file_dir=self.lineCalibrantDir.text(),
            calibrant_file_formulae=self.checkCalibrantFormulae.isChecked(),
            calibrant_search_window=self.spinSearchWindow.value(),
            calibrant_min_intsy=self.spinMinIntsy.value(),
        )

    def get_input_files(
        self,
        validate_paths: bool = False,
    ) -> list[Path | str]:
        """
        Retrieves input filepaths from UI, and optionally
        ensures that they exist/valid
        :return:
        """

        current_filepaths: list[str] = [
            self.listInputFiles.item(x).text() for x in
            range(self.listInputFiles.count())
        ]

        if not validate_paths:
            return current_filepaths

    def execute(self):
        """
        Executes the main conversion script in a separate thread
        (after checking that some parameters have been specified)
        :return:
        """
        self.retrieve_parameters_from_ui()

        if not self.parameters:
            return None

        # Create worker + thread and connect signals
        self.thread = QThread()
        self.worker = ProcessingWorker(
            params=self.parameters,
        )
        self.worker.moveToThread(self.thread)


        self.thread.started.connect(self.worker.run)

        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.finished.connect(self.process_ended)

        self.worker.progress.connect(self.update_progress)

        self.worker.error.connect(self.error)

        # Start thread
        self.thread.start()

        # Disable execute btn
        self.btnExecute.setEnabled(False)
        self.thread.finished.connect(
            lambda: self.btnExecute.setEnabled(True)
        )

        # Create progress dialog
        self.progress_dialog = QtWidgets.QProgressDialog(
            'Processing ...',
            "Cancel",
            0,
            100,
            self,
        )
        self.progress_dialog.setWindowTitle(
            'Processing .mzML files'
        )
        self.progress_dialog.setModal(True)
        self.progress_dialog.show()

    def update_progress(
        self,
        value: int,
    ):
        if self.progress_dialog:
            self.progress_dialog.setValue(value)

    def process_ended(
        self,
    ):
        if self.progress_dialog:
            self.progress_dialog.close()
            self.progress_dialog.deleteLater()
            self.progress_dialog = None

    def error(
        self,
        msg: str,
    ):
        QtWidgets.QErrorMessage(self).showMessage(msg)
        self.process_ended()


    def btnReportFilepathPressed(self):
        """
        Called when user clicks to load a report .json file
        """
        report_filepath, _ = QtWidgets.QFileDialog.getOpenFileName(
            parent=self,
            caption='Open report .json file (generated after calibrating)',
            filter="Calibration report file (*.json)"
        )

        if not report_filepath:
            return

        self.lineReportFilepath.setText(report_filepath)

    def btnReportPlotPressed(self):

        report_filepath = self.lineReportFilepath.text()

        if (not Path(report_filepath).exists() or
            Path(report_filepath).suffix.lower() != '.json'):
            return

        print('Displaying plot..')
        win: pg.GraphicsLayoutWidget = calibration.read_and_display_calibration_report(
            report_path=Path(report_filepath),
            return_graphics_layout_widget=True,
        )

        # Delete old plot
        # See: https://stackoverflow.com/questions/4528347/
        for i in range(self.layoutCalibrationPlot.count()):
            _ = self.layoutCalibrationPlot.itemAt(i).widget()
            if _:
                self.layoutCalibrationPlot.removeWidget(_)

        # Add new plot
        self.layoutCalibrationPlot.addWidget(win)
