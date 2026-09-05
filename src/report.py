"""Phase 4: the figures and the machine-readable report.

Figure 1: calibration curves, Condition A (solid) vs B (dashed), CI bands, diagonal,
abstention threshold line, MAF-stratified small multiples.
Figure 2: abstention trade-off - threshold sweep vs retained accuracy, per condition.
README numbers are generated from benchmark_report.json, never hand-typed.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pydantic import BaseModel

from src.abstention import AbstentionOutcome
from src.calibration import CalibrationResult, bootstrap_mean_ci


class BinReport(BaseModel):
    lo: float
    hi: float
    n: int
    mean_dr2: float
    mean_r2: float
    r2_ci_low: float
    r2_ci_high: float
    gap: float


class ConditionReport(BaseModel):
    condition: str
    n_variants: int
    n_zero_variance_excluded: int
    n_missing_in_imputed: int
    ece: float
    bins: list[BinReport]
    maf_strata: dict[str, list[BinReport]]
    maf_strata_ece: dict[str, float]
    abstention: dict


class BenchmarkReport(BaseModel):
    seed: int
    n_bins: int
    bootstrap_reps: int
    accuracy_floor: float
    threshold_selected_on: str
    conditions: list[ConditionReport]


def _bins_to_report(c: CalibrationResult) -> list[BinReport]:
    return [BinReport(**vars(b)) for b in c.bins]


def condition_report(
    condition: str,
    calib: CalibrationResult,
    strata: dict[str, CalibrationResult],
    abst: AbstentionOutcome,
    n_zero_var: int,
    n_missing: int,
) -> ConditionReport:
    return ConditionReport(
        condition=condition,
        n_variants=calib.n_variants,
        n_zero_variance_excluded=n_zero_var,
        n_missing_in_imputed=n_missing,
        ece=round(calib.ece, 6),
        bins=_bins_to_report(calib),
        maf_strata={k: _bins_to_report(v) for k, v in strata.items()},
        maf_strata_ece={k: round(v.ece, 6) for k, v in strata.items()},
        abstention=vars(abst),
    )


def _plot_curve(ax, calib: CalibrationResult, label: str, style: str):
    xs = [b.mean_dr2 for b in calib.bins]
    ys = [b.mean_r2 for b in calib.bins]
    lo = [b.r2_ci_low for b in calib.bins]
    hi = [b.r2_ci_high for b in calib.bins]
    ax.plot(xs, ys, style, label=label, marker="o", markersize=3)
    ax.fill_between(xs, lo, hi, alpha=0.2)


def figure_calibration(
    calib_a: CalibrationResult,
    calib_b: CalibrationResult,
    strata_a: dict[str, CalibrationResult],
    strata_b: dict[str, CalibrationResult],
    threshold: float | None,
    out_path: Path,
):
    keys = sorted(set(strata_a) | set(strata_b))
    fig, axes = plt.subplots(1, 1 + len(keys), figsize=(4.2 * (1 + len(keys)), 4.2))
    axes = np.atleast_1d(axes)

    panels = [("All evaluation variants", calib_a, calib_b)] + [
        (k.replace("maf_", "MAF ").replace("_", " to "), strata_a.get(k), strata_b.get(k))
        for k in keys
    ]
    for ax, (title, ca, cb) in zip(axes, panels, strict=True):
        ax.plot([0, 1], [0, 1], color="grey", lw=0.8, ls=":")
        if ca:
            _plot_curve(ax, ca, "A: matched (AFR panel)", "-")
        if cb:
            _plot_curve(ax, cb, "B: mismatched (EUR panel)", "--")
        if threshold is not None:
            ax.axvline(threshold, color="black", lw=0.8, ls="-.",
                       label=f"abstain < {threshold:.1f}")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("mean DR2 (promised)")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    axes[0].set_ylabel("mean empirical r$^2$ (observed)")
    axes[0].legend(fontsize=7, loc="upper left")
    fig.suptitle("Is DR2 calibrated, and does calibration survive panel mismatch?", fontsize=11)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def figure_tradeoff(
    dr2_a: np.ndarray, r2_a: np.ndarray,
    dr2_b: np.ndarray, r2_b: np.ndarray,
    bootstrap_reps: int, seed: int, out_path: Path,
):
    thresholds = np.linspace(0.0, 0.95, 20)
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    for name, dr2, r2, style in (
        ("A: matched", dr2_a, r2_a, "-"),
        ("B: mismatched", dr2_b, r2_b, "--"),
    ):
        fr_abst, acc, lo, hi = [], [], [], []
        for t in thresholds:
            mask = dr2 >= t
            if mask.sum() < 10:
                break
            fr_abst.append(1 - mask.mean())
            vals = r2[mask]
            acc.append(vals.mean())
            ci_lo, ci_hi = bootstrap_mean_ci(vals, bootstrap_reps, seed=seed)
            lo.append(ci_lo)
            hi.append(ci_hi)
        ax.plot(fr_abst, acc, style, label=name, marker="o", markersize=3)
        ax.fill_between(fr_abst, lo, hi, alpha=0.2)
    ax.set_xlabel("fraction of variants abstained")
    ax.set_ylabel("mean empirical r$^2$ among retained")
    ax.set_title("Abstention trade-off per condition", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
