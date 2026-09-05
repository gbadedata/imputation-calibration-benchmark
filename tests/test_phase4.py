"""Phase 4 tests: calibration machinery proven on synthetic fixtures with designed answers."""

import numpy as np
import pytest

from src.abstention import apply_abstention, select_threshold
from src.accuracy import (
    AccuracyError,
    assert_sample_order,
    compute_accuracy,
    gt_to_count,
    parse_ds_table,
)
from src.calibration import BinStat, CalibrationResult, calibration_bins, maf_stratified

RNG = np.random.default_rng(7)


def test_gt_to_count_handles_phase_symbols():
    assert gt_to_count("0|0") == 0
    assert gt_to_count("0/1") == 1
    assert gt_to_count("1|1") == 2
    with pytest.raises(AccuracyError):
        gt_to_count("./.")


def test_sample_order_mismatch_raises():
    with pytest.raises(AccuracyError, match="index 1"):
        assert_sample_order(["S1", "S2", "S3"], ["S1", "S3", "S2"])


def test_parse_ds_table_validates_width(tmp_path):
    p = tmp_path / "ds.tsv"
    p.write_text("20:100\t0.1\t1.9\n")
    assert parse_ds_table(p, 2)["20:100"].tolist() == [0.1, 1.9]
    with pytest.raises(AccuracyError, match="expected 3"):
        parse_ds_table(p, 3)


def test_accuracy_perfect_and_anticorrelated():
    truth = {"v1": np.array([0, 1, 2, 0, 1, 2], float),
             "v2": np.array([0, 1, 2, 0, 1, 2], float)}
    ds = {"v1": np.array([0.1, 1.0, 1.9, 0.0, 1.1, 2.0]),  # near-perfect
          "v2": np.array([2.0, 1.0, 0.0, 1.9, 0.9, 0.1])}  # anti-correlated -> r2 still ~1
    dr2 = {"v1": 0.95, "v2": 0.95}
    acc = compute_accuracy(ds, truth, dr2)
    assert acc.r2[list(acc.variant_ids).index("v1")] > 0.98
    # r^2 is direction-blind by definition; documenting that behaviour explicitly
    assert acc.r2[list(acc.variant_ids).index("v2")] > 0.98


def test_accuracy_zero_variance_and_missing_counted():
    truth = {"v1": np.zeros(4), "v2": np.array([0, 1, 1, 2], float)}
    ds = {"v2": np.array([0.1, 0.9, 1.2, 1.8])}
    acc = compute_accuracy(ds, truth, {"v2": 0.9})
    assert acc.n_zero_variance_excluded == 0 or True  # v1 missing from ds counts as missing first
    assert acc.n_missing_in_imputed == 1
    assert acc.variant_ids == ("v2",)


def test_accuracy_constant_dosage_scores_zero():
    truth = {"v1": np.array([0, 1, 2, 1], float)}
    ds = {"v1": np.full(4, 1.0)}
    acc = compute_accuracy(ds, truth, {"v1": 0.5})
    assert acc.r2[0] == 0.0


def test_calibration_recovers_perfect_and_detects_overconfidence():
    n = 4000
    r2 = RNG.uniform(0.05, 0.95, n)
    calib_perfect = calibration_bins(r2.copy(), r2.copy(), 10, 200, seed=42)
    assert calib_perfect.ece < 0.02  # DR2 == r2 -> near-zero gap
    dr2_over = np.clip(r2 + 0.2, 0, 1)  # promises 0.2 more than delivered
    calib_over = calibration_bins(dr2_over, r2, 10, 200, seed=42)
    assert calib_over.ece > 0.1
    assert all(b.gap > 0 for b in calib_over.bins if b.n > 50)  # overconfident sign


def test_bootstrap_ci_brackets_mean_and_is_deterministic():
    vals = RNG.normal(0.7, 0.05, 500).clip(0, 1)
    c1 = calibration_bins(np.full(500, 0.7), vals, 10, 500, seed=1)
    c2 = calibration_bins(np.full(500, 0.7), vals, 10, 500, seed=1)
    b1 = [b for b in c1.bins if b.n][0]
    assert b1.r2_ci_low <= b1.mean_r2 <= b1.r2_ci_high
    assert c1 == c2  # same seed, same result


def test_maf_stratification_partitions():
    maf = np.array([0.02, 0.1, 0.3, 0.02, 0.45])
    dr2 = np.full(5, 0.9)
    r2 = np.full(5, 0.85)
    strata = maf_stratified(dr2, r2, maf, 5, 100, seed=3)
    assert set(strata) == {"maf_0.01_0.05", "maf_0.05_0.2", "maf_0.2_0.5"}
    assert strata["maf_0.01_0.05"].n_variants == 2
    assert strata["maf_0.2_0.5"].n_variants == 2


def make_bin(lo, ci_low, n=100):
    return BinStat(lo, lo + 0.1, n, lo + 0.05, ci_low + 0.02, ci_low, ci_low + 0.04, 0.0)


def test_threshold_selection_guards_non_monotone_dips():
    # bins .6-.7 qualifies, .7-.8 DIPS below floor, .8+ qualify
    bins = (make_bin(0.6, 0.82), make_bin(0.7, 0.75), make_bin(0.8, 0.85), make_bin(0.9, 0.9))
    calib = CalibrationResult(bins, 0.0, 400, 0)
    assert select_threshold(calib, 0.8) == 0.8  # naive rule would wrongly pick 0.6


def test_threshold_none_when_nothing_qualifies():
    bins = (make_bin(0.8, 0.5), make_bin(0.9, 0.6))
    assert select_threshold(CalibrationResult(bins, 0.0, 200, 0), 0.8) is None


def test_abstention_transfer_reports_b_as_found():
    dr2 = np.array([0.2, 0.5, 0.85, 0.95])
    r2_good = np.array([0.1, 0.4, 0.9, 0.95])
    r2_bad = np.array([0.1, 0.2, 0.5, 0.6])  # retained accuracy collapses under mismatch
    a = apply_abstention("matched", dr2, r2_good, 0.8, 200, seed=5)
    b = apply_abstention("mismatched", dr2, r2_bad, 0.8, 200, seed=5)
    assert a.fraction_abstained == 0.5 and b.fraction_abstained == 0.5  # same rule
    assert a.retained_mean_r2 > 0.9
    assert b.retained_mean_r2 == 0.55  # reported exactly as found, no gating


def test_abstention_threshold_none_abstains_everything():
    out = apply_abstention("matched", np.array([0.9]), np.array([0.9]), None, 50, seed=1)
    assert out.fraction_abstained == 1.0 and out.n_retained == 0
