"""
Functions for reading/writing data to/from disk
"""
from pathlib import Path

import pyopenms as oms
from molmass import Formula


def read_mzml(
    input_path: Path,
) -> oms.MSExperiment:
    """
    Reads .mzML files as MSExperiment() objects
    """
    input_exp = oms.MSExperiment()
    oms.MzMLFile().load(
        filename=str(input_path),
        in_1=input_exp,
    )
    return input_exp


def load_calibrant_series(
    file_path: Path,
    formulae: bool = False,
) -> list[float]:
    """
    Reads the file in `file_path` and extracts a list of m/z values
    that can be used for calibration. File must be one m/z or formula per line
    Parameters
    ----------
    file_path:  Path to calibrant series file (one m/z or formula per line)
    formulae:   Whether the calibrant file is a list of formulae
                    (will calculate exact masses if True; make sure to
                    include charge, i.e. C6H12O6+)

    Returns
    -------
    List of m/z values to be used for calibration
    """
    if not file_path.exists():
        raise FileNotFoundError(
            f"File not found at {file_path}"
        )

    calibrant_masses: list[float] = []
    with open(file_path) as calibrant_file:
        for line in calibrant_file:
            if not line:
                # Skip empty lines
                continue

            match formulae:
                case False:
                    calibrant_masses.append(
                        float(line)
                    )
                case True:
                    calibrant_masses.append(
                        Formula(line).monoisotopic_mass
                    )

    return sorted(calibrant_masses)
