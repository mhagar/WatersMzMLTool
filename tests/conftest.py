import pytest
import pyopenms as oms
from src.file_utils import read_mzml

from pathlib import Path

# Path to test file
TEST_FILE_PATH = Path(
    "/home/mh/PycharmProjects/WatersMzMLTool/tests/20251008_QC_9.mzML"
)

@pytest.fixture
def test_file_path() -> Path:
    return TEST_FILE_PATH


@pytest.fixture
def ms_exp(test_file_path) -> oms.MSExperiment:
    exp = read_mzml(test_file_path)
    return exp
