"""
Scripts specifically for handling quirks in .MzML files from Waters
"""
import re
from copy import deepcopy
from typing import Literal

import pyopenms as oms


PATTERN: re.Pattern = re.compile(
    r'^function=(\d+)'
)

def fix_missing_ms_level_labels(
    ms_experiment: oms.MSExperiment,
    experiment_type: Literal['dda', 'dia'],
) -> oms.MSExperiment:
    """
    Given an MSExperiment, returns a copy where the MSLevel
    attributes are adjusted according to the 'Native ID'

    Parameters
    ----------
    ms_experiment:      oms.MSExperiment
    experiment_type:    Either 'dda' or 'dia'. Note: dda not supported yet

    Returns
    -------
    an MSExperiment where the spectra have adjusted MSLevels
    """
    adjusted_exp: oms.MSExperiment = deepcopy(ms_experiment)
    adjusted_spectra: list[oms.MSSpectrum] = []
    for spectrum in ms_experiment.getSpectra():
        # Get spectrum info
        native_id: str = spectrum.getNativeID()
        ms_level: int = native_id_to_ms_level(
            native_id=native_id,
            experiment_type=experiment_type,
        )

        spectrum.setMSLevel(ms_level)

        adjusted_spectra.append(
            spectrum
        )

    adjusted_exp.setSpectra(adjusted_spectra)

    return adjusted_exp

def apply_centroiding(
    ms_experiment: oms.MSExperiment,
) -> oms.MSExperiment:
    """
    Wrapper around OpenMS's PeakPickerIterative.
    # TODO: Add parameters!

    Parameters
    ----------
    ms_experiment:      oms.MSExperiment


    Returns
    -------
    an MSExperiment containing centroided spectra
    """
    adjusted_exp: oms.MSExperiment = deepcopy(ms_experiment)
    # adjusted_spectra: list[oms.MSSpectrum] = []

    picker_iterative = oms.PeakPickerIterative()

    picker_iterative.pickExperiment(
        input=ms_experiment,
        output=adjusted_exp,
    )

    return adjusted_exp


def native_id_to_ms_level(
    native_id: str,
    experiment_type: Literal['dia', 'dda'],
) -> int:
    """
    Given a native_id of the format: 'function=x process=y scan=z',
    returns x as an integer.

    Parameters
    ----------
    native_id:      String of the format 'function=x process=y scan=z'
    experiment_type:Either 'dia' or 'dda'. Note: only 'dia' supported currently

    Returns
    -------
    Integer
    """
    # Get function number
    hit = PATTERN.search(
        native_id
    )

    if not hit:
        raise ValueError(
            f"Unable to parse native id: {native_id}"
        )

    function_num: int = int(hit.group(1))

    match experiment_type:
        case 'dda':
            raise NotImplementedError(
                "Haven't implemented DDA yet"
            )

        case 'dia':
            return function_num
