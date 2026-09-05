"""Command-line entry points used by scripts/phase1.sh.

Usage:
    python -m src.cli design-samples --panel <panel_file> --outdir <dir> [--removals <file>]
    python -m src.cli make-sites --mode thin --site-table <tsv> --outdir <dir>
    python -m src.cli qc-report --qc-log <log> --snplist <f> --king-id <f> --outdir <dir>
    python -m src.cli bgen-verify --original <sig> --round-tripped <sig>

All thresholds and seeds come from config.settings.Settings (ICB_* env overridable);
the CLI deliberately exposes no tuning flags so runs are reproducible from the repo alone.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from config.settings import Settings
from src.bgen_roundtrip import read_signature, verify_round_trip
from src.imputation_check import check_condition, parse_dr2_table, read_lines
from src.phase4 import run_phase4
from src.qc_report import build_report, render_markdown
from src.sample_design import design_samples, read_panel, read_sample_list, write_design
from src.sites import (
    partition_sites,
    partition_sites_by_thinning,
    read_hm3_rsids,
    read_site_table,
    write_partition,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="icb")
    sub = parser.add_subparsers(dest="command", required=True)

    p_design = sub.add_parser("design-samples")
    p_design.add_argument("--panel", type=Path, required=True)
    p_design.add_argument("--outdir", type=Path, required=True)
    p_design.add_argument("--removals", type=Path, default=None)

    p_sites = sub.add_parser("make-sites")
    p_sites.add_argument("--mode", choices=["thin", "hm3"], required=True)
    p_sites.add_argument("--site-table", type=Path, required=True)
    p_sites.add_argument("--outdir", type=Path, required=True)
    p_sites.add_argument("--hm3", type=Path, default=None, help="required in hm3 mode")

    p_qc = sub.add_parser("qc-report")
    p_qc.add_argument("--qc-log", type=Path, required=True)
    p_qc.add_argument("--snplist", type=Path, required=True)
    p_qc.add_argument("--king-id", type=Path, required=True)
    p_qc.add_argument("--outdir", type=Path, required=True)

    p_bgen = sub.add_parser("bgen-verify")
    p_bgen.add_argument("--original", type=Path, required=True)
    p_bgen.add_argument("--round-tripped", type=Path, required=True)

    p_imp = sub.add_parser("impute-check")
    p_imp.add_argument("--condition", choices=["matched", "mismatched"], required=True)
    p_imp.add_argument("--header", type=Path, required=True)
    p_imp.add_argument("--samples", type=Path, required=True)
    p_imp.add_argument("--targets", type=Path, required=True)
    p_imp.add_argument("--marker-ids", type=Path, required=True)
    p_imp.add_argument("--eval-ids", type=Path, required=True)
    p_imp.add_argument("--dr2-table", type=Path, required=True)
    p_imp.add_argument("--out", type=Path, required=True)

    p_cal = sub.add_parser("calibrate")
    p_cal.add_argument("--extract-dir", type=Path, required=True,
                       help="dir with {A,B}_ds.tsv, {A,B}_dr2.tsv, truth_gt.tsv, *_samples.txt")
    p_cal.add_argument("--report", type=Path, required=True)
    p_cal.add_argument("--figures-dir", type=Path, required=True)

    args = parser.parse_args(argv)
    settings = Settings()

    if args.command == "design-samples":
        panel = read_panel(args.panel)
        removals = read_sample_list(args.removals) if args.removals else set()
        design = design_samples(panel, removals, settings)
        summary = write_design(design, args.outdir, settings)
    elif args.command == "make-sites":
        records = read_site_table(args.site_table)
        if args.mode == "hm3":
            if args.hm3 is None:
                parser.error("--hm3 is required in hm3 mode")
            hm3 = read_hm3_rsids(args.hm3)
            partition = partition_sites(records, hm3, settings.eval_maf_floor)
        else:
            partition = partition_sites_by_thinning(
                records, settings.n_array_sites, settings.eval_maf_floor
            )
        summary = write_partition(partition, args.outdir)
    elif args.command == "qc-report":
        report = build_report(
            args.qc_log.read_text(),
            snplist_path=args.snplist,
            king_id_path=args.king_id,
            geno=settings.geno_missing_max,
            maf=settings.maf_min,
            hwe=settings.hwe_p,
            king_cutoff=settings.king_cutoff,
        )
        args.outdir.mkdir(parents=True, exist_ok=True)
        (args.outdir / "qc_report.json").write_text(report.model_dump_json(indent=2) + "\n")
        (args.outdir / "qc_report.md").write_text(render_markdown(report) + "\n")
        summary = report.model_dump()
    elif args.command == "bgen-verify":
        verdict = verify_round_trip(
            read_signature(args.original), read_signature(args.round_tripped)
        )
        summary = {"verdict": verdict}
    elif args.command == "calibrate":
        summary = run_phase4(args.extract_dir, args.report, args.figures_dir, settings)
    else:  # impute-check
        dr2_values, n_missing = parse_dr2_table(args.dr2_table)
        metrics = check_condition(
            condition=args.condition,
            header_text=args.header.read_text(),
            sample_ids=read_lines(args.samples),
            target_ids=read_lines(args.targets),
            output_marker_ids=set(read_lines(args.marker_ids)),
            eval_site_ids=set(read_lines(args.eval_ids)),
            dr2_values=dr2_values,
            n_dr2_missing=n_missing,
            settings=settings,
        )
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(metrics.model_dump_json(indent=2) + "\n")
        summary = metrics.model_dump()

    json.dump(summary, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
