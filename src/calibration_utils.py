"""
Functions for calibrating Waters .MzML files
"""
from typing import Iterable, Literal, Optional
from copy import deepcopy

import numpy as np
import pyopenms as oms
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

from src.ms_utils import get_spectrum_type, find_signal


def accumulate_corrected_spectra(
    spectra: list[oms.MSSpectrum],
    fit_func: callable,
    corr_type: Literal["absolute", "relative"],
    calibrant_series: list[float],
    min_intsy: float,
    search_window: float,
    include_intensities: bool,
    inplace: bool = False,
    min_num_cal_signals: int = 3,  # TODO: Expose to user
) -> list[oms.MSSpectrum]:
    """
	Iterates through an MSExperiment's spectra and looks for lockmass spectra.
	Uses the calibrants defined by `calibrant_series` to correct the
	acquisitions between lockmass scans (as well as the lockmass scans
	themselves).

	Parameters
	----------
	spectra:            List of oms.MSSpectrum objects
	fit_func:           Mathematical function to fit to mass error
	corr_type:          Whether the function should work on ppm or Da error
	calibrant_series:   Series of theoretical m/z values to calibrate with
	search_window:      m/z window when finding calibrant signals
	min_num_cal_signals:Minimum number of calibrant signals that must be
	                    detected in a calibration spectrum for it to be
	                    valid. Default: 3
	min_intsy:          Minimum signal intensity when finding calibrant signals

	include_intensities: Whether `fit_func` expects just m/z values, or
						 also intensity values
	inplace:            Whether to adjust the input objects directly, return a
						copy

	Returns
	-------
	A list of oms.MSSpectrum objects with corrected m/z values
	"""
    if not inplace:
        spectra = deepcopy(spectra)

    corrected_spectra: list[oms.MSSpectrum] = []
    calibration_spectra: list[oms.MSSpectrum] = []
    most_recent_cal_params: Iterable = []

    for spectrum in spectra:
        # If calibration spectrum, update calibration params
        if get_spectrum_type(spectrum) == "cal":
            cal_params, num_cal_signals = get_calibration_params(
                cal_spec=spectrum,
                fit_func=fit_func,
                corr_type=corr_type,
                calibrant_series=calibrant_series,
                search_window=search_window,
                min_intsy=min_intsy,
                include_intensities=include_intensities,
            )

            if num_cal_signals > min_num_cal_signals:
                most_recent_cal_params = cal_params
                calibration_spectra.append(spectrum)

        if len(most_recent_cal_params) == 0:
            # No lockmass scan encountered yet.
            # Just add the spectrum to the stack uncorrected and move on
            corrected_spectra.append(spectrum)
            continue

        # Format spectrum object into spectrum array
        mz_arr, intsy_arr = spectrum.get_peaks()
        spec_arr = np.zeros((len(mz_arr), 2))

        spec_arr[:, 0] = mz_arr
        spec_arr[:, 1] = intsy_arr

        # Apply calibration
        calibrated_spec_arr = apply_calibration(
            spec=spec_arr,
            fit_func=fit_func,
            params=most_recent_cal_params,
            corr_type=corr_type,
            include_intensities=include_intensities,
        )

        # Fix spectrum
        spectrum.set_peaks(
            (calibrated_spec_arr[:, 0], calibrated_spec_arr[:, 1])
        )

        corrected_spectra.append(spectrum)

    return corrected_spectra


def apply_global_spline_correction(
    spectra: list[oms.MSSpectrum],
    corr_type: Literal["absolute", "relative"],
    calibrant_series: list[float],
    search_window: float,
    min_intsy: float,
    smoothing_factor: float = 0.0,
    spline_k: int = 3,  # Default cubic spline
    return_spline: bool = False,
    inplace: bool = False,
) -> tuple[list[oms.MSSpectrum], Optional[UnivariateSpline]]:
    """
    Applies a correction where *all* calibrant scans are considered.
    This is intended to be used as a second pass, after doing
    chunk-by-chunk calibration.
    """
    if not inplace:
        spectra = deepcopy(spectra)

    # Get all the calibration spectra
    calibration_spectra: list[oms.MSSpectrum] = []
    for spectrum in spectra:
        if get_spectrum_type(spectrum) == "cal":
            calibration_spectra.append(spectrum)

    # Retrieve average calibrant signals in those arrays
    calmzs = []
    for calmass in calibrant_series:
        this_calibrant = []  # To be averaged
        for spectrum in calibration_spectra:
            hit_mz, hit_intsy = find_signal(
                target_mz=calmass,
                spec=spectrum,
                window_da=search_window,
                min_intsy=min_intsy,
            )

            if not hit_mz:
                continue

            match corr_type:
                case "absolute":
                    error = calmass - hit_mz
                case "relative":
                    error = (calmass - hit_mz) / (calmass * 1e-6)

            this_calibrant.append((hit_mz, error, hit_intsy))

        if not this_calibrant:
            continue

        this_calibrant = np.array(this_calibrant)
        calmzs.append(
            (
                this_calibrant[:, 0].mean(),
                this_calibrant[:, 1].mean(),
                this_calibrant[:, 2].mean(),
            )
        )


    # Sort then turn into array
    calmzs = sorted(
        calmzs,
        key=lambda x: x[0]
    )
    calmzs = np.array(calmzs)

    # Get calibration params using combined signals
    calibration_spline = get_spline(
        mz_values=calmzs[:, 0],
        error_values=calmzs[:, 1],
        smoothing_factor=smoothing_factor,
        spline_k=spline_k,
    )

    # Now apply calibration to all the spectra
    corrected_spectra: list[oms.MSSpectrum] = []
    for spectrum in spectra:
        # Format spectrum object into spectrum array
        mz_arr, intsy_arr = spectrum.get_peaks()
        spec_arr = np.zeros((len(mz_arr), 2))
        spec_arr[:, 0] = mz_arr
        spec_arr[:, 1] = intsy_arr

        corrected_spec_arr = apply_spline_calibration(
            spec=spec_arr,
            spline=calibration_spline,
            corr_type=corr_type,
        )

        # Fix spectrum
        spectrum.set_peaks(
            (
                corrected_spec_arr[:, 0],
                corrected_spec_arr[:, 1],
            )
        )

        corrected_spectra.append(spectrum)

    if return_spline:
        return corrected_spectra, calibration_spline

    return corrected_spectra, None


def get_calibration_params(
    cal_spec: oms.MSSpectrum,
    fit_func: callable,
    corr_type: Literal["absolute", "relative"],
    calibrant_series: list[float],
    search_window: float,
    min_intsy: float,
    include_intensities: bool = True,
) -> tuple[Iterable, int]:
    """
    Given an MSSpectrum object, finds calibrant signals and performs a fit
    using `fit_func`. If `include_intensities` is True, also fits the
    signal intensities

    Returns a tuple containing:
        curve fitting parameters
        number of calibrant signals used
    """
    calmzs: list = []
    for calmass in calibrant_series:
        hit_mz, hit_intsy = find_signal(
            target_mz=calmass,
            spec=cal_spec,
            window_da=search_window,
            min_intsy=min_intsy,
        )

        if not hit_mz:
            continue

        match corr_type:
            case "absolute":
                error = calmass - hit_mz
            case "relative":
                error = (calmass - hit_mz) / (calmass * 1e-6)

        calmzs.append((hit_mz, error, hit_intsy))

    calmzs: np.ndarray = np.array(sorted(calmzs))
    if calmzs.size == 0:
        return [], calmzs.shape[0]

    if include_intensities:
        params = curve_fit(
                fit_func,
                calmzs[:, [0, 2]],  # mz, intsy
                calmzs[:, 1],  # error (mass discrepancy)
            )[0]
    else:
        params = curve_fit(
            fit_func,
            calmzs[:, 0],  # mz only
            calmzs[:, 1],  # error (mass discrepnacy)
        )[0]


    return params, calmzs.shape[0]


def apply_calibration(
    spec: np.ndarray,
    fit_func: callable,
    params: Iterable,
    corr_type: Literal["absolute", "relative"],
    include_intensities: bool,
) -> np.ndarray:

    if len(list(params)) == 0: raise ValueError("No params given")

    mz = spec[:, 0]

    if include_intensities:
        dependent_variable = spec  # mz, intsy
    else:
        dependent_variable = mz

    corr_factor = fit_func(
        dependent_variable,
        *params,
    )

    match corr_type:
        case "absolute":
            corr_mz: np.ndarray = mz + corr_factor
        case "relative":
            corr_mz: np.ndarray = mz + (1e-6 * mz * (corr_factor - 100))

    result = spec.copy()
    result[:, 0] = corr_mz

    return result


def get_spline(
    mz_values: np.ndarray,
    error_values: np.ndarray,
    smoothing_factor: float = 0.0,
    spline_k: int = 3,  # Default cubic spline
) -> UnivariateSpline:
    """
    Given arrays of mz_values and error_values, returns a UnivariateSpline
    object that can be used to interpolate

    Parameters:
    -----------
    mz_values : np.ndarray
        Array of m/z values to fit spline with

    error_values: np.ndarray
        Array of error values to fit spline with

    smoothing_factor : float
        Smoothing factor for the spline (0 = exact interpolation)
    spline_k : int
        Degree of the spline (1=linear, 2=quadratic, 3=cubic)

    Returns:
    --------
    UnivariateSpline

    """
    return UnivariateSpline(
        x=mz_values,
        y=error_values,
        k=spline_k,
        s=smoothing_factor
    )


def apply_spline_calibration(
    spec: np.ndarray,
    spline: UnivariateSpline,
    corr_type: Literal["absolute", "relative"],
) -> np.ndarray:
    """
    Adjusts a spec array using a UnivariateSpline
    (prepared using get_spline() ),
    """
    mz = spec[:, 0]
    corr_factor = spline(
        mz
    )

    match corr_type:
        case "absolute":
            corr_mz: np.ndarray = mz + corr_factor
        case "relative":
            corr_mz: np.ndarray = mz + (1e-6 * mz * (corr_factor - 100))

    result = spec.copy()
    result[:, 0] = corr_mz

    return result


def remove_lockmass_scans(
    spectra: list[oms.MSSpectrum]
) -> list[oms.MSSpectrum]:
    return [
        spec for spec in spectra
        if get_spectrum_type(spec) != 'cal'
    ]