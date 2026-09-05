"""Phase 3 validation: check each Beagle output before it is allowed into Phase 4.

Inputs are plain-text extracts produced by scripts/phase3.sh with bcftools (header,
sample list, marker IDs, evaluation-site DR2 table); this module never parses VCF
directly, keeping it unit-testable on synthetic fixtures.

Gating is asymmetric by design (BUILD-SPEC Phase 3): the matched condition (A) must
meet strict floors; the mismatched condition (B) has its metrics REPORTED, with only a
catastrophic coverage floor enforced - dropped AFR-common sites under a EUR panel are a
mismatch phenomenon being measured, not a pipeline failure to gate away.
"""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel

from config.settings import Settings


class ImputationCheckError(ValueError):
    """Raised when a Beagle output fails validation for its condition."""


class ImputationMetrics(BaseModel):
    condition: str  # "matched" | "mismatched"
    dr2_declared: bool
    ds_declared: bool
    n_samples: int
    samples_match_targets: bool
    n_markers_output: int
    n_eval_sites_expected: int
    n_eval_sites_covered: int
    eval_coverage: float
    n_dr2_missing: int
    dr2_missing_fraction: float
    n_dr2_out_of_range: int


def header_declares(header_text: str, field_id: str) -> bool:
    return f"<ID={field_id}," in header_text


def read_lines(path: Path) -> list[str]:
    return [ln.strip() for ln in Path(path).read_text().splitlines() if ln.strip()]


def parse_dr2_table(path: Path) -> tuple[dict[str, float], int]:
    """ID<TAB>DR2 rows -> ({id: dr2}, n_missing). '.' or empty DR2 counts as missing.
    Number=A note: biallelic-only design means one value; a comma-separated value is
    treated as malformed for this pipeline and raises."""
    values: dict[str, float] = {}
    n_missing = 0
    for i, ln in enumerate(Path(path).read_text().splitlines(), start=1):
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) != 2:
            raise ImputationCheckError(f"dr2 table line {i}: expected 2 fields")
        vid, raw = parts
        if raw in (".", ""):
            n_missing += 1
            continue
        if "," in raw:
            raise ImputationCheckError(
                f"dr2 table line {i}: multi-allele DR2 '{raw}' - biallelic-only design violated"
            )
        values[vid] = float(raw)
    return values, n_missing


def check_condition(
    *,
    condition: str,
    header_text: str,
    sample_ids: list[str],
    target_ids: list[str],
    output_marker_ids: set[str],
    eval_site_ids: set[str],
    dr2_values: dict[str, float],
    n_dr2_missing: int,
    settings: Settings,
) -> ImputationMetrics:
    if condition not in ("matched", "mismatched"):
        raise ImputationCheckError(f"unknown condition {condition!r}")

    dr2_ok = header_declares(header_text, "DR2")
    ds_ok = header_declares(header_text, "DS")
    samples_match = sorted(sample_ids) == sorted(target_ids)

    covered = eval_site_ids & output_marker_ids
    coverage = len(covered) / len(eval_site_ids) if eval_site_ids else 0.0

    n_eval_with_dr2 = len(dr2_values) + n_dr2_missing
    missing_frac = n_dr2_missing / n_eval_with_dr2 if n_eval_with_dr2 else 1.0
    out_of_range = sum(1 for v in dr2_values.values() if not (0.0 <= v <= 1.0))

    m = ImputationMetrics(
        condition=condition,
        dr2_declared=dr2_ok,
        ds_declared=ds_ok,
        n_samples=len(sample_ids),
        samples_match_targets=samples_match,
        n_markers_output=len(output_marker_ids),
        n_eval_sites_expected=len(eval_site_ids),
        n_eval_sites_covered=len(covered),
        eval_coverage=round(coverage, 6),
        n_dr2_missing=n_dr2_missing,
        dr2_missing_fraction=round(missing_frac, 6),
        n_dr2_out_of_range=out_of_range,
    )

    # Universal requirements (both conditions)
    problems = []
    if not dr2_ok:
        problems.append("DR2 not declared in header")
    if not ds_ok:
        problems.append("DS not declared in header")
    if not samples_match:
        problems.append("output samples do not match target list")
    if out_of_range:
        problems.append(f"{out_of_range} DR2 values outside [0, 1]")

    # Condition-specific floors
    if condition == "matched":
        if coverage < settings.coverage_floor_matched:
            problems.append(
                f"matched coverage {coverage:.4f} below floor {settings.coverage_floor_matched}"
            )
        if missing_frac > settings.dr2_missing_max_matched:
            problems.append(
                f"matched DR2-missing fraction {missing_frac:.4f} above "
                f"{settings.dr2_missing_max_matched}"
            )
    else:
        if coverage < settings.coverage_floor_mismatched:
            problems.append(
                f"mismatched coverage {coverage:.4f} below catastrophic floor "
                f"{settings.coverage_floor_mismatched}"
            )

    if problems:
        raise ImputationCheckError(
            f"condition {condition} failed validation: " + "; ".join(problems)
        )
    return m
