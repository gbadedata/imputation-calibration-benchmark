"""Phase 4: abstention analysis.

Threshold rule (BUILD-SPEC 4.3, stated choice not clinical standard): choose the
smallest DR2 bin lower-edge such that EVERY bin at or above it has an observed-accuracy
bootstrap CI lower bound >= the accuracy floor (0.8 by default). The "every bin above"
form guards against non-monotone dips that the naive smallest-qualifying-bin rule would
step over. Threshold is chosen on Condition A ONLY, then applied unchanged to B; the
outcome on B is reported exactly as found.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src.calibration import CalibrationResult, bootstrap_mean_ci


@dataclass(frozen=True)
class AbstentionOutcome:
    condition: str
    threshold: float | None  # None = no qualifying bin; abstain on everything
    fraction_abstained: float
    n_retained: int
    retained_mean_r2: float | None
    retained_r2_ci_low: float | None
    retained_r2_ci_high: float | None


def select_threshold(calib_a: CalibrationResult, accuracy_floor: float) -> float | None:
    """Smallest bin lower-edge with all bins at/above clearing the CI-lower floor."""
    bins = calib_a.bins
    for i, b in enumerate(bins):
        if all(bb.r2_ci_low >= accuracy_floor for bb in bins[i:]):
            return b.lo
    return None


def apply_abstention(
    condition: str,
    dr2: np.ndarray,
    r2: np.ndarray,
    threshold: float | None,
    bootstrap_reps: int,
    seed: int,
) -> AbstentionOutcome:
    n = len(dr2)
    if threshold is None:
        return AbstentionOutcome(condition, None, 1.0, 0, None, None, None)
    mask = dr2 >= threshold
    n_ret = int(mask.sum())
    if n_ret == 0:
        return AbstentionOutcome(condition, threshold, 1.0, 0, None, None, None)
    retained = r2[mask]
    ci_lo, ci_hi = bootstrap_mean_ci(retained, bootstrap_reps, seed=seed)
    return AbstentionOutcome(
        condition=condition,
        threshold=threshold,
        fraction_abstained=round(1.0 - n_ret / n, 6),
        n_retained=n_ret,
        retained_mean_r2=round(float(retained.mean()), 6),
        retained_r2_ci_low=round(ci_lo, 6),
        retained_r2_ci_high=round(ci_hi, 6),
    )
