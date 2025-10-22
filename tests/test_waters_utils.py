from src.waters_utils import fix_missing_ms_level_labels, apply_centroiding
import pytest

def test_fix_missing_ms_level_labels(ms_exp):
    fixed_exp = fix_missing_ms_level_labels(
        ms_experiment=ms_exp,
        experiment_type='dia',
    )

    # Accumulate scan ms-levels
    ms_levels = []
    for spec in fixed_exp.getSpectra():
        ms_levels.append(
            spec.getMSLevel()
        )

    assert len(set(ms_levels)) > 1, (
        f'After fixing MS level labels, experiment had these ms levels: '
        f'{set(ms_levels)}'
    )


def test_apply_centroiding(ms_exp):
    centroided_exp = apply_centroiding(
        ms_experiment=ms_exp
    )

    assert ms_exp.size() == centroided_exp.size(), (
        "Experiments have unequal scan numbers after centroiding"
    )

    # Get the middle spectrum
    mid_scan = centroided_exp.size() // 2

    centroided_exp_test_scan = centroided_exp.getSpectra()[mid_scan]
    profile_exp_test_scan = ms_exp.getSpectra()[mid_scan]

    assert centroided_exp_test_scan.size() < profile_exp_test_scan.size(), (
        "Centroided scan should have *less* signals than profile scan"
    )


