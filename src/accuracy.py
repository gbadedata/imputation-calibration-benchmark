"""Phase 4: per-variant empirical imputation accuracy against withheld truth.

Inputs are bcftools-query extracts (never raw VCF):
- DS table:    ID<TAB>ds_1<TAB>...<TAB>ds_100   (imputed ALT dosage per target)
- truth table: ID<TAB>gt_1<TAB>...<TAB>gt_100   (GT strings from truth_targets)
- two sample-order files, asserted identical before any correlation is computed

Per variant: empirical accuracy = squared Pearson correlation between imputed dosage
and true genotype (0/1/2) across the 100 targets. Variants with zero truth variance in
the target cohort are excluded and counted (BUILD-SPEC 4.1) - never silently dropped.
Truth MAF (target-cohort, not global AF) is computed here for stratification.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


class AccuracyError(ValueError):
    """Raised on malformed or misaligned Phase 4 inputs."""


@dataclass(frozen=True)
class AccuracyResult:
    variant_ids: tuple[str, ...]  # variants with computed accuracy, aligned arrays below
    r2: np.ndarray  # empirical r^2 per variant
    dr2: np.ndarray  # Beagle DR2 per variant
    maf: np.ndarray  # truth target-cohort MAF per variant
    n_zero_variance_excluded: int
    n_missing_in_imputed: int  # eval sites absent from the imputed output (e.g. the dup)


def assert_sample_order(truth_samples: list[str], imputed_samples: list[str]) -> None:
    if truth_samples != imputed_samples:
        pairs = zip(truth_samples, imputed_samples, strict=False)
        first = next(i for i, (x, y) in enumerate(pairs) if x != y)
        raise AccuracyError(
            "sample order differs between truth and imputed extracts; correlations "
            f"would be scrambled. First divergence at index {first}"
        )


def parse_ds_table(path: Path, n_samples: int) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for i, ln in enumerate(Path(path).read_text().splitlines(), start=1):
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) != n_samples + 1:
            raise AccuracyError(f"DS table line {i}: {len(parts)-1} values, expected {n_samples}")
        out[parts[0]] = np.array([float(x) for x in parts[1:]], dtype=np.float64)
    return out


def gt_to_count(gt: str) -> int:
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or any(a not in ("0", "1") for a in alleles):
        raise AccuracyError(f"unexpected GT {gt!r} (biallelic non-missing expected)")
    return alleles.count("1")


def parse_truth_table(path: Path, n_samples: int) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for i, ln in enumerate(Path(path).read_text().splitlines(), start=1):
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) != n_samples + 1:
            raise AccuracyError(
                f"truth table line {i}: {len(parts)-1} values, expected {n_samples}"
            )
        out[parts[0]] = np.array([gt_to_count(g) for g in parts[1:]], dtype=np.float64)
    return out


def parse_dr2_map(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    for ln in Path(path).read_text().splitlines():
        if not ln.strip():
            continue
        vid, raw = ln.split("\t")
        if raw not in (".", ""):
            out[vid] = float(raw)
    return out


def compute_accuracy(
    ds: dict[str, np.ndarray],
    truth: dict[str, np.ndarray],
    dr2: dict[str, float],
) -> AccuracyResult:
    eval_ids = sorted(truth.keys())
    n_missing_imputed = sum(1 for v in eval_ids if v not in ds or v not in dr2)

    ids: list[str] = []
    r2s: list[float] = []
    dr2s: list[float] = []
    mafs: list[float] = []
    n_zero_var = 0

    for vid in eval_ids:
        if vid not in ds or vid not in dr2:
            continue
        t = truth[vid]
        d = ds[vid]
        if np.var(t) == 0.0:
            n_zero_var += 1
            continue
        p = float(np.mean(t)) / 2.0
        maf = min(p, 1.0 - p)
        if np.var(d) == 0.0:
            # constant dosage against varying truth: correlation undefined, accuracy 0
            r2 = 0.0
        else:
            r = float(np.corrcoef(d, t)[0, 1])
            r2 = r * r
        ids.append(vid)
        r2s.append(r2)
        dr2s.append(dr2[vid])
        mafs.append(maf)

    return AccuracyResult(
        variant_ids=tuple(ids),
        r2=np.array(r2s),
        dr2=np.array(dr2s),
        maf=np.array(mafs),
        n_zero_variance_excluded=n_zero_var,
        n_missing_in_imputed=n_missing_imputed,
    )
