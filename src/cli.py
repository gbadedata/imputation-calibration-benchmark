"""Command-line entry points used by scripts/phase1.sh.

Usage:
    python -m src.cli design-samples --panel <panel_file> --outdir <dir> [--removals <file>]
    python -m src.cli make-sites --hm3 <w_hm3.snplist> --site-table <tsv> --outdir <dir>

All thresholds and seeds come from config.settings.Settings (ICB_* env overridable);
the CLI deliberately exposes no tuning flags so runs are reproducible from the repo alone.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from config.settings import Settings
from src.sample_design import design_samples, read_panel, read_sample_list, write_design
from src.sites import partition_sites, read_hm3_rsids, read_site_table, write_partition


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="icb")
    sub = parser.add_subparsers(dest="command", required=True)

    p_design = sub.add_parser("design-samples")
    p_design.add_argument("--panel", type=Path, required=True)
    p_design.add_argument("--outdir", type=Path, required=True)
    p_design.add_argument("--removals", type=Path, default=None)

    p_sites = sub.add_parser("make-sites")
    p_sites.add_argument("--hm3", type=Path, required=True)
    p_sites.add_argument("--site-table", type=Path, required=True)
    p_sites.add_argument("--outdir", type=Path, required=True)

    args = parser.parse_args(argv)
    settings = Settings()

    if args.command == "design-samples":
        panel = read_panel(args.panel)
        removals = read_sample_list(args.removals) if args.removals else set()
        design = design_samples(panel, removals, settings)
        summary = write_design(design, args.outdir, settings)
    else:
        hm3 = read_hm3_rsids(args.hm3)
        records = read_site_table(args.site_table)
        partition = partition_sites(records, hm3, settings.eval_maf_floor)
        summary = write_partition(partition, args.outdir)

    json.dump(summary, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
