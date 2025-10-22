"""
Functions for parsing MSExperiment and MSSpectrum objects
"""
from typing import Literal

import numpy as np
import pyopenms as oms


def get_spectrum_type(
    spec: oms.MSSpectrum,
) -> Literal["ms1", "ms2", "cal", None]:
    """
    TODO: This assumes DIA (i.e. three scan functions)
    Parameters
    ----------
    spec:   MSSpectrum object

    Returns
    -------
    Literal["ms1", "ms2", "cal"]
    """
    spec_id = spec.getNativeID()

    if "function=1" in spec_id:
        return "ms1"

    if "function=2" in spec_id:
        # return "cal"
        return "ms2"

    if "function=3" in spec_id:
        return "cal"

    if not spec_id:
        return None

    raise ValueError(
        f"Unable to determine spectrum type given spec_id: {spec_id}"
    )


def find_signal(
    target_mz: float,
    spec: oms.MSSpectrum,
    window_da: float,
    min_intsy: float = 0.0,
) -> tuple[float, float] | tuple[None, None]:
    """
    Given a spectrum array, returns tuple (mass, intensity) corresponding to
    the tallest signal within (window_ppm) of (target_mz), OR an empty
    array if no match.

    Optionally, min_intsy describes threshold at which a signal counts
    as being 'detected'.
    """
    idx: int = spec.findHighestInWindow(
        mz=target_mz,
        tolerance_left=window_da,
        tolerance_right=window_da,
    )

    if idx == -1:
        # None match
        return None, None

    mz, intsy = spec.get_mz_array()[idx], spec.get_intensity_array()[idx]

    if intsy < min_intsy:
        # None match
        return None, None

    return mz, intsy


def calc_rsme(
    errors: np.ndarray,
) -> float:
    return np.sqrt(
        np.sum(
            errors ** 2
        ) / (
            errors.size
        )
    )


def find_mass_errors(
    spectra: list[oms.MSSpectrum],
    theoretical_mzs: list[float],
    window_da: float,
    min_intsy: float,
    spectrum_type: Literal['ms1', 'ms2', 'cal'],
) -> list[list[tuple[float, float]]]:
    """
	Given a list of MSSpectrum objects and a list of theoretical mz values,
	searches for the requested signals (within the constraints defined by
	min_intsy, spectrum_type, and window_da)

	Parameters
	----------
	window_da
	spectra
	theoretical_mzs
	min_intsy
	spectrum_type

	Returns
	-------

	"""
    errors_vs_mz: list[list[tuple[float, float]]] = []
    for spectrum in spectra:
        spectrum: oms.MSSpectrum
        if get_spectrum_type(spectrum) != spectrum_type: continue

        # Accumulate mass error for all calibrants in this scan
        err_vs_mz: list[tuple[float, float]] = []
        for theoretical in theoretical_mzs:
            observed_mz, intsy = find_signal(
                target_mz=theoretical,
                spec=spectrum,
                min_intsy=min_intsy,
                window_da=window_da,
            )

            if not observed_mz: continue

            err_vs_mz.append(
                (
                    observed_mz - theoretical,
                    theoretical,
                )
            )

        errors_vs_mz.append(
            err_vs_mz
        )

    return errors_vs_mz


