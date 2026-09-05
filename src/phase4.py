"""Phase 4 orchestrator: extracts -> accuracy -> calibration -> abstention -> outputs.

Expected files in extract_dir (produced by scripts/phase4.sh):
  truth_gt.tsv, truth_samples.txt, A_ds.tsv, B_ds.tsv, A_dr2.tsv, B_dr2.tsv,
  A_samples.txt, B_samples.txt
"""

from __future__ import annotations

import json
from pathlib import Path

from config.settings import Settings
from src.abstention import apply_abstention, select_threshold
from src.accuracy import (
    assert_sample_order,
    compute_accuracy,
    parse_dr2_map,
    parse_ds_table,
    parse_truth_table,
)
from src.calibration import calibration_bins, maf_stratified
from src.report import (
    BenchmarkReport,
    condition_report,
    figure_calibration,
    figure_tradeoff,
)


def run_phase4(
    extract_dir: Path, report_path: Path, figures_dir: Path, settings: Settings
) -> dict:
    d = Path(extract_dir)
    truth_samples = [s for s in (d / "truth_samples.txt").read_text().split() if s]
    n = len(truth_samples)
    truth = parse_truth_table(d / "truth_gt.tsv", n)

    results = {}
    for cond_id, cond_name in (("A", "matched"), ("B", "mismatched")):
        imp_samples = [s for s in (d / f"{cond_id}_samples.txt").read_text().split() if s]
        assert_sample_order(truth_samples, imp_samples)
        ds = parse_ds_table(d / f"{cond_id}_ds.tsv", n)
        dr2 = parse_dr2_map(d / f"{cond_id}_dr2.tsv")
        acc = compute_accuracy(ds, truth, dr2)
        calib = calibration_bins(
            acc.dr2, acc.r2, settings.dr2_bins, settings.bootstrap_reps, settings.seed
        )
        strata = maf_stratified(
            acc.dr2, acc.r2, acc.maf, settings.dr2_bins,
            settings.bootstrap_reps, settings.seed,
        )
        results[cond_id] = dict(name=cond_name, acc=acc, calib=calib, strata=strata)

    threshold = select_threshold(
        results["A"]["calib"], settings.abstention_accuracy_floor
    )
    condition_reports = []
    for cond_id in ("A", "B"):
        r = results[cond_id]
        abst = apply_abstention(
            r["name"], r["acc"].dr2, r["acc"].r2, threshold,
            settings.bootstrap_reps, settings.seed,
        )
        condition_reports.append(
            condition_report(
                r["name"], r["calib"], r["strata"], abst,
                r["acc"].n_zero_variance_excluded, r["acc"].n_missing_in_imputed,
            )
        )

    report = BenchmarkReport(
        seed=settings.seed,
        n_bins=settings.dr2_bins,
        bootstrap_reps=settings.bootstrap_reps,
        accuracy_floor=settings.abstention_accuracy_floor,
        threshold_selected_on="matched (Condition A) only; transferred unchanged to B",
        conditions=condition_reports,
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report.model_dump_json(indent=2) + "\n")

    figures_dir = Path(figures_dir)
    figure_calibration(
        results["A"]["calib"], results["B"]["calib"],
        results["A"]["strata"], results["B"]["strata"],
        threshold, figures_dir / "fig1_calibration.png",
    )
    figure_tradeoff(
        results["A"]["acc"].dr2, results["A"]["acc"].r2,
        results["B"]["acc"].dr2, results["B"]["acc"].r2,
        settings.bootstrap_reps, settings.seed, figures_dir / "fig2_tradeoff.png",
    )

    return {
        "threshold": threshold,
        "A": {"ece": report.conditions[0].ece, "n": report.conditions[0].n_variants,
              "abstained": report.conditions[0].abstention["fraction_abstained"],
              "retained_r2": report.conditions[0].abstention["retained_mean_r2"]},
        "B": {"ece": report.conditions[1].ece, "n": report.conditions[1].n_variants,
              "abstained": report.conditions[1].abstention["fraction_abstained"],
              "retained_r2": report.conditions[1].abstention["retained_mean_r2"]},
        "report": str(report_path),
        "figures": [str(figures_dir / "fig1_calibration.png"),
                    str(figures_dir / "fig2_tradeoff.png")],
    }


def load_report(path: Path) -> dict:
    return json.loads(Path(path).read_text())
