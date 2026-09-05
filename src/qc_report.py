"""Phase 2 QC report: parse plink2 outputs into a validated, consistency-checked report.

Parsing philosophy (BUILD-SPEC 2.3): plink2 log wording can drift between builds, so the
parser is tolerant (a missing removal line reads as 0) but the report generator enforces a
HARD accounting identity before writing anything:

    n_input_variants - (geno + maf + hwe removals) == n_kept_variants (from the snplist)

If the identity fails - because a wording changed and a removal went uncounted, or because
plink2 applied a filter we did not model - the report raises instead of shipping wrong
numbers. Silent zeroes cannot survive this check.
"""

from __future__ import annotations

import re
from pathlib import Path

from pydantic import BaseModel, model_validator


class QCReportError(ValueError):
    """Raised when plink2 outputs cannot be reconciled into a consistent report."""


_PATTERNS = {
    "n_input_variants": [
        re.compile(r"(\d+) variants loaded from"),
        re.compile(r"--vcf: (\d+) variants scanned"),
    ],
    "n_samples": [re.compile(r"(\d+) samples \(")],
    "n_geno_removed": [re.compile(r"(\d+) variants? removed due to missing genotype data")],
    "n_maf_removed": [
        re.compile(r"(\d+) variants? removed due to allele frequency threshold"),
        re.compile(r"(\d+) variants? removed due to minor allele threshold"),
    ],
    "n_hwe_removed": [re.compile(r"(\d+) variants? removed due to Hardy-Weinberg exact test")],
}


def parse_plink_log(text: str) -> dict[str, int]:
    """Extract counts from a plink2 .log. Removal lines default to 0 when absent;
    the accounting identity in QCReport catches any wording drift that hides a removal."""
    out: dict[str, int] = {}
    for key, patterns in _PATTERNS.items():
        value: int | None = None
        for pat in patterns:
            m = pat.search(text)
            if m:
                value = int(m.group(1))
                break
        if value is None:
            if key in ("n_input_variants", "n_samples"):
                raise QCReportError(
                    f"could not find {key} in plink2 log; inspect the log and update parser"
                )
            value = 0
        out[key] = value
    return out


def count_id_file(path: Path) -> int:
    """Count sample IDs in a plink2 .id file (e.g. king.cutoff.out.id), skipping '#' headers."""
    return sum(
        1
        for ln in Path(path).read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    )


def count_snplist(path: Path) -> int:
    return sum(1 for ln in Path(path).read_text().splitlines() if ln.strip())


class QCReport(BaseModel):
    n_input_variants: int
    n_samples: int
    n_geno_removed: int
    n_maf_removed: int
    n_hwe_removed: int
    n_kept_variants: int
    n_kinship_removed_samples: int
    geno_threshold: float
    maf_threshold: float
    hwe_threshold: float
    king_cutoff: float

    @model_validator(mode="after")
    def accounting_identity(self) -> QCReport:
        removed = self.n_geno_removed + self.n_maf_removed + self.n_hwe_removed
        expected_kept = self.n_input_variants - removed
        if expected_kept != self.n_kept_variants:
            raise QCReportError(
                f"QC accounting failed: {self.n_input_variants} input - {removed} removed "
                f"= {expected_kept}, but snplist has {self.n_kept_variants}. A filter went "
                "uncounted (log wording drift?) or an unmodelled filter ran. Inspect the log."
            )
        return self


def build_report(
    log_text: str,
    snplist_path: Path,
    king_id_path: Path,
    *,
    geno: float,
    maf: float,
    hwe: float,
    king_cutoff: float,
) -> QCReport:
    counts = parse_plink_log(log_text)
    return QCReport(
        **counts,
        n_kept_variants=count_snplist(snplist_path),
        n_kinship_removed_samples=count_id_file(king_id_path),
        geno_threshold=geno,
        maf_threshold=maf,
        hwe_threshold=hwe,
        king_cutoff=king_cutoff,
    )


def render_markdown(r: QCReport) -> str:
    return "\n".join(
        [
            "# Phase 2 QC report",
            "",
            f"Kinship: {r.n_kinship_removed_samples} samples removed at KING cutoff "
            f"{r.king_cutoff} (2nd degree), across all AFR+EUR candidates, before the redraw.",
            "",
            "Chip-cohort variant QC (targets only - see BUILD-SPEC 2.2 for why panels are "
            "not used in chip QC):",
            "",
            "| Step | Threshold | Variants removed |",
            "|---|---|---|",
            f"| Input | - | {r.n_input_variants} |",
            f"| Call rate (--geno) | {r.geno_threshold} | {r.n_geno_removed} |",
            f"| MAF (--maf) | {r.maf_threshold} | {r.n_maf_removed} |",
            f"| HWE (--hwe) | {r.hwe_threshold} | {r.n_hwe_removed} |",
            f"| **Kept** | - | **{r.n_kept_variants}** |",
            "",
            f"Samples in QC run: {r.n_samples}.",
            "",
            "Accounting identity (input - removals == kept) verified at report build time.",
        ]
    )
