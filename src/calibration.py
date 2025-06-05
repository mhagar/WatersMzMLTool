"""
Functions for applying calibrations to an oms.MSExperiment object
"""
import pyopenms as oms
import pyqtgraph as pg
import numpy as np
from PyQt5.QtCore import Qt
from scipy.interpolate import UnivariateSpline

import src.calibration_utils as cal_utils
import src.calibration_funcs as cal_funcs
import src.ms_utils as ms_utils

import json
from pathlib import Path
from typing import Optional


def correct_ms_exp_first_pass(
    spectra: list[oms.MSSpectrum],
    calibrant_series: list[float],
    calibrant_window_da: float,
    min_intsy: float,
) -> list[oms.MSSpectrum]:
    """
    Given an MSExperiment object, searches the lockmass scans for the masses
    given by `calibrant_series`, fits a calibration curve (parabola),
    then corrects all the scans accordingly.

    This is the first pass - a second pass should be used as well.

    Parameters
    ----------
    exp 					MSExperiment object
    calibrant_series 		List of calibrant m/z values to search for
    min_intsy 				Minimum signal intensity to consider
    calibrant_window_da 	Search window

    Returns
    -------
    An adjusted MSExperiment
    """
    # Apply first pass correction (i.e. parabola function)
    corrected_spectra = cal_utils.accumulate_corrected_spectra(
        spectra=spectra,
        fit_func=cal_funcs.parabola,
        # fit_func=cal_funcs.parabola_w_intensity_term,
        corr_type="absolute",
        calibrant_series=calibrant_series,
        search_window=calibrant_window_da,
        min_intsy=min_intsy,
        include_intensities=False,
        # include_intensities=True,
    )

    return corrected_spectra


def correct_ms_exp_second_pass(
    spectra: list[oms.MSSpectrum],
    calibrant_series: list[float],
    calibrant_window_da: float,
    min_intsy: float,
) -> tuple[list[oms.MSSpectrum], UnivariateSpline]:
    """
    Given an MSExperiment object, searches the lockmass scans for the masses
    given by `calibrant_series`, fits a spline-based calibration curve,
    then corrects all the scans accordingly.

    This is the second pass - should be used after the parabola-based
    calibration.

    Parameters
    ----------
    exp 					MSExperiment object
    calibrant_series 		List of calibrant m/z values to search for
    min_intsy 				Minimum signal intensity to consider
    calibrant_window_da 	Search window

    Returns
    -------
    An adjusted MSExperiment
    """
    corrected_spectra, spline = cal_utils.apply_global_spline_correction(
        spectra=spectra,
        corr_type='absolute',
        calibrant_series=calibrant_series,
        search_window=calibrant_window_da,
        min_intsy=min_intsy,
        return_spline=True,
    )
    return corrected_spectra, spline


def generate_calibration_report(
    uncorrected_spectra: list[oms.MSSpectrum],
    corrected_spectra_first_pass: list[oms.MSSpectrum],
    corrected_spectra_secnd_pass: list[oms.MSSpectrum],
    spline: UnivariateSpline,
    calibrants: list[float],
    window_da: float,
    min_intsy: float,
) -> dict:
    """
    Generates a report that details the calibrant signals used
    (uncorrected, after first pass correction, and after second pass)

    Outputs a dictionary that can be written as a JSON file

    The format looks like:
    >>> {
    >>> 	"calibration_scan_1": {         # Calibration scan number
    >>> 		"uncorrected": [
    >>> 			[556.54, 0.023],
    >>> 			[785.53, 0.054]
    >>>
    >>> 		],
    >>> 		"first_corr": [
    >>> 			[556.54, 0.004],
    >>> 			[785.54, 0.003],
    >>> 		],
    >>> 		"second_corr": [
    >>> 			[556.54, 0.002],
    >>> 			[785.54, 0.001],
    >>> 		],
    >>> 	},
    >>> 	"calibration_scan_2": {
    >>> 	# ... (A leaf for each calibration scan)
    >>> 	},
    >>>     "parameters" : {
    >>>         "search_window": 0.1,
    >>>         "minimum_intensity": 1000,
    >>>         "first_corr_params": {
    >>>         },
    >>>         "second_corr_params": {
    >>>         }
    >>>     }
    >>> }

    :param calibrants: A list of theoretical calibrant m/zs
    :param window_da: The m/z tolerance to search for them with
    :param min_intsy: The minimum intensity a calibrant signal must be
    :param uncorrected_spectra: A list of MSSpectrum objects pre-correction
    :param corrected_spectra_first_pass: A list of MSSpectrum objs
    :param corrected_spectra_secnd_pass: A list of MSSpectrum objs
    :param spline: The UnivariateSpline object used to perform second calibtn
    :return: A dictionary that can be written as a JSON file
    """

    output = {
        'calibration_scans': {},
        'parameters': {},
    }

    mass_errors = {
        'uncorrected': [],
        'first_corr': [],
        'second_corr': [],
    }

    for label, spectra in zip(
        ('uncorrected', 'first_corr', 'second_corr'),
        (uncorrected_spectra, corrected_spectra_first_pass, corrected_spectra_secnd_pass,)
    ):
        mass_errors[label] = ms_utils.find_mass_errors(
            spectra=spectra,
            spectrum_type='cal',
            theoretical_mzs=calibrants,
            window_da=window_da,
            min_intsy=min_intsy,
        )

    # Add calibration_scan_N data
    for scan_num in range(len(mass_errors['uncorrected'])):
        cal_scan_num = f"cal_scan_{scan_num}"
        output['calibration_scans'][cal_scan_num] = {}

        for label in ['uncorrected', 'first_corr', 'second_corr']:

            output['calibration_scans'][cal_scan_num][label]: list[list[float, float]] =[
                [
                    float(x[0]),
                    float(x[1]),
                ] for x in mass_errors[label][scan_num]
            ]

    # Now add parameters data
    output['parameters'] = {
        'calibrants': [float(x) for x in calibrants],
        'search_window': window_da,
        'minimum_intensity': min_intsy,
        'first_corr_params': {},
        'second_corr_params': [float(x) for x in spline(calibrants)],
    }

    return output


def read_and_display_calibration_report(
    report_path: Path,
    return_graphics_layout_widget: bool = False,
):
    """
    Given a path to a report .JSON file, displays a calibration error plot
    :param report_path: Path to .json file
    :param return_graphics_layout_widget: Whether to display window, or return widget
    :return:
    """
    with open(report_path, 'r') as f:
        report: dict = json.load(f)

    (uncorrected_mass_errors,
    first_corr_mass_errors,
    second_corr_mass_errors) = [], [], []

    for mass_errors, label in zip(
        (uncorrected_mass_errors,
         first_corr_mass_errors,
         second_corr_mass_errors),
        ('uncorrected',
         'first_corr',
         'second_corr')
    ):
        for scan_num, data in report['calibration_scans'].items():
            arr = np.array(data[label])
            if arr.any():
                mass_errors.append(arr)

    spline = np.array(report['parameters']['second_corr_params'])

    return show_calibrant_error_plot(
        uncorrected_mass_errors=uncorrected_mass_errors,
        first_corr_mass_errors=first_corr_mass_errors,
        second_corr_mass_errors=second_corr_mass_errors,
        spline=spline,
        return_graphics_layout_widget=return_graphics_layout_widget,
    )


def show_calibrant_error_plot(
    uncorrected_mass_errors: list[np.ndarray],
    first_corr_mass_errors: list[np.ndarray],
    second_corr_mass_errors: list[np.ndarray],
    spline: np.ndarray,
    return_graphics_layout_widget: bool = False,
) -> Optional['pg.GraphicsLayoutWidget']:
    """
    Given a list of MSSpectrum objects, searches lockmass scans
    for m/z values given by `calibrants`, and plots a 'm/z vs mass error'
    graph
    """
    # Initialize graphicslayoutwidget
    win = pg.GraphicsLayoutWidget(
        show=True
    )

    uncorrected_plot: pg.PlotWidget = win.addPlot(row=0, col=0)
    first_pass_plot: pg.PlotWidget = win.addPlot(row=1, col=0)
    second_pass_plot: pg.PlotWidget = win.addPlot(row=2, col=0)

    first_pass_plot.setXLink(uncorrected_plot)
    second_pass_plot.setXLink(uncorrected_plot)

    first_pass_plot.setYLink(second_pass_plot)


    # Populate the plots
    for pw, mz_vs_errors, label in zip(
        [uncorrected_plot, first_pass_plot, second_pass_plot],
        [uncorrected_mass_errors, first_corr_mass_errors, second_corr_mass_errors],
        ["Uncorrected data", "First-pass correction", "Second-pass correction"]
    ):
        mz_vs_errors: list[np.ndarray]
        populate_plot_widget(
            pw=pw,
            mz_vs_errors=mz_vs_errors,
        )

        rmse = ms_utils.calc_rsme(
            errors=np.array(mz_vs_errors)[:, :, 0].flatten()
        )

        pw.setTitle(
            f"{label}: {len(mz_vs_errors)} "
            f"calibration scans : "
            f"RMSE = {rmse:.5f} Da"
        )

        # Add spline to first pass correction plot
        if label == "First-pass correction":
            mz_values = np.array(mz_vs_errors[0])[:, 1]

            pw.addItem(
                pg.PlotDataItem(
                    mz_values,
                    -spline,
                    pen=pg.mkPen(
                        (255, 0, 255, 200),
                        width=2,
                    ),
                    name="Spline used for second pass"
                )
            )

    if return_graphics_layout_widget:
        return win

    else:
        pg.exec()


def populate_plot_widget(
    pw: pg.PlotWidget,
    mz_vs_errors: list[np.ndarray],
):
    pw.addLegend(
        offset=(10, 10)
    )

    # Add the actual error plots:
    for idx, mz_vs_error in enumerate(mz_vs_errors):
        if not mz_vs_error.any():
            continue

        pw.addItem(
            pg.PlotDataItem(
                mz_vs_error[:, 1],
                mz_vs_error[:, 0],
                pen=pg.intColor(idx, alpha=60),
                symbolPen=pg.intColor(idx, alpha=60),
                symbolBrush=pg.intColor(idx, alpha=60),
                symbol='x',
            )
        )

    # Add the '1 ppm' indicator
    mz_values = np.array(mz_vs_errors[-1])[:, 1],
    mz_range = np.linspace(
        start=np.min(mz_values),
        stop=np.max(mz_values),
        num=100,
    )

    ppm_error = mz_range * 1e-6
    for sign in (-1, +1):
        pw.addItem(
            pg.PlotDataItem(
                mz_range,
                ppm_error * sign,
                pen=pg.mkPen(
                    (0, 255, 255, 128),
                    width=1,
                    style=Qt.DashLine,
                ),
                name=f"{sign * 1} ppm"
            )
        )

    # Add a 0 line
    pw.addItem(
        pg.InfiniteLine(
            pos=0,
            angle=0,
            pen=pg.mkPen(
                (128, 128, 128)
            )
        )
    )

    # Set axis labels
    pw.setLabel(
        "left",
        "Mass error (Da)"
    )

    pw.setLabel(
        "bottom",
        "m/z",
    )