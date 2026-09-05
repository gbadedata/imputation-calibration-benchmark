"""Phase 4: DR2 calibration analysis with bootstrap confidence intervals.

Question: within each DR2 bin, does mean DR2 (Beagle's promise) match mean empirical
r^2 (what the truth shows)? Every reported bin statistic carries a percentile bootstrap
CI over variants (BUILD-SPEC 4.2: no headline number ships without an interval).

ECE definition (falsifiable, per module docstring requirement): the variant-count-
weighted mean of |mean(DR2) - mean(r^2)| across non-empty bins.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MAF_STRATA = ((0.01, 0.05), (0.05, 0.20), (0.20, 0.50))


@dataclass(frozen=True)
class BinStat:
    lo: float
    hi: float
    n: int
    mean_dr2: float
    mean_r2: float
    r2_ci_low: float
    r2_ci_high: float
    gap: float  # mean_dr2 - mean_r2 (positive = overconfident)


@dataclass(frozen=True)
class CalibrationResult:
    bins: tuple[BinStat, ...]
    ece: float
    n_variants: int
    n_empty_bins: int


def bootstrap_mean_ci(
    values: np.ndarray, reps: int, seed: int, alpha: float = 0.05
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(values)
    idx = rng.integers(0, n, size=(reps, n))
    means = values[idx].mean(axis=1)
    return float(np.quantile(means, alpha / 2)), float(np.quantile(means, 1 - alpha / 2))


def calibration_bins(
    dr2: np.ndarray,
    r2: np.ndarray,
    n_bins: int,
    bootstrap_reps: int,
    seed: int,
) -> CalibrationResult:
    if len(dr2) != len(r2):
        raise ValueError("dr2 and r2 length mismatch")
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    bins: list[BinStat] = []
    n_empty = 0
    total_gap_weight = 0.0
    for k in range(n_bins):
        lo, hi = float(edges[k]), float(edges[k + 1])
        mask = (dr2 >= lo) & (dr2 < hi) if k < n_bins - 1 else (dr2 >= lo) & (dr2 <= hi)
        n = int(mask.sum())
        if n == 0:
            n_empty += 1
            continue
        sub_r2 = r2[mask]
        m_dr2 = float(dr2[mask].mean())
        m_r2 = float(sub_r2.mean())
        ci_lo, ci_hi = bootstrap_mean_ci(sub_r2, bootstrap_reps, seed=seed + k)
        gap = m_dr2 - m_r2
        total_gap_weight += abs(gap) * n
        bins.append(BinStat(lo, hi, n, m_dr2, m_r2, ci_lo, ci_hi, gap))
    n_total = int(sum(b.n for b in bins))
    ece = total_gap_weight / n_total if n_total else 0.0
    return CalibrationResult(tuple(bins), float(ece), n_total, n_empty)


def maf_stratified(
    dr2: np.ndarray,
    r2: np.ndarray,
    maf: np.ndarray,
    n_bins: int,
    bootstrap_reps: int,
    seed: int,
) -> dict[str, CalibrationResult]:
    out: dict[str, CalibrationResult] = {}
    for lo, hi in MAF_STRATA:
        mask = (maf >= lo) & (maf < hi) if hi < 0.5 else (maf >= lo) & (maf <= hi)
        if mask.sum() == 0:
            continue
        key = f"maf_{lo}_{hi}"
        out[key] = calibration_bins(dr2[mask], r2[mask], n_bins, bootstrap_reps, seed)
    return out
