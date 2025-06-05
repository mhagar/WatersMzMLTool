"""
A Worker thread for running the script in GUI mode,
and giving progress updates
"""
from waters_mzml_tools import process_mzml_file_cli, _setup_parser

from PyQt5.QtCore import QObject, pyqtSignal
from pyqtgraph import GraphicsLayoutWidget

from pathlib import Path
import traceback
import argparse
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from gui_mode import Parameters


class ProcessingWorker(QObject):
    # Signals to communicate w main thread
    progress = pyqtSignal(int)
    finished = pyqtSignal()
    error = pyqtSignal(str)

    def __init__(
        self,
        params: 'Parameters',
    ):
        super().__init__()
        self.params = params

    def run(self):
        """
        Runs in background thread
        :return:
        """
        try:
            self.process_mzml_files(self.params)
            self.finished.emit()

        except Exception as e:
            full_traceback = traceback.format_exc()
            self.error.emit(
                full_traceback
                # str(e)
            )

    def process_mzml_files(
        self,
        params: 'Parameters',
    ):
        """
        Main entrypoint into mzml processing script
        :param params:
        :return:
        """
        # Assemble arguments except for --input:
        parser: argparse.ArgumentParser = _setup_parser()
        args = parser.parse_args(
            params.format_as_cli_arguments()
        )

        total_num_files: int = len(params.input_files)

        for idx, filepath in enumerate(params.input_files):
            print(
                f"Processing: {filepath}"
            )
            process_mzml_file_cli(
                input_path=Path(filepath),
                args=args,
            )

            self.progress.emit(
	            int(100 * idx/total_num_files)
            )



