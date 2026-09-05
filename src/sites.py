"""Site design: partition chr20 biallelic SNPs into array sites and evaluation sites.

Two modes:

hm3 mode (retained, unused for the v5b data): intersect rsIDs with the HapMap3 list.
  The 1kGP 20130502 v5b chr20 release carries NO rsIDs (ID column is '.' on every row,
  verified 2026-09-04), so this mode is infeasible for the primary data. The code and
  tests remain because the policy is correct and the reason it is unused is documented.

thin mode (primary): deterministic uniform thinning to n_array sites. IDs are chr:pos
  (set upstream by bcftools annotate --set-id). The thinning pool is restricted to
  common variants (maf_floor <= AF <= 1 - maf_floor) BEFORE selection: real genotyping
  arrays are deliberately common-biased, and an "array" of rare singletons would give
  the imputation engine nothing to phase on. Selection is evenly spaced by position
  order - reproducible from the repo alone, no seed involved.

Shared policy, in order (each exclusion counted, never silent - BUILD-SPEC Phase 1.2):
1. rows with missing IDs ('.' or empty) are excluded from BOTH site sets
2. duplicate IDs keep the first occurrence; later occurrences excluded from both
3. array/evaluation split per mode; array and evaluation are disjoint by construction
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path


class SiteTableError(ValueError):
    """Raised on malformed site-table or hm3 input."""


@dataclass(frozen=True)
class SiteRecord:
    rsid: str
    pos: int
    af: float


@dataclass(frozen=True)
class SitePartition:
    array_sites: tuple[SiteRecord, ...]
    evaluation_sites: tuple[SiteRecord, ...]
    n_missing_id: int
    n_duplicate_id: int
    n_eval_excluded_maf: int


def read_hm3_rsids(path: Path) -> frozenset[str]:
    lines = Path(path).read_text().splitlines()
    if not lines:
        raise SiteTableError(f"hm3 snplist is empty: {path}")
    header = lines[0].split()
    if not header or header[0] != "SNP":
        raise SiteTableError(f"unexpected hm3 header {lines[0]!r}; expected to start with 'SNP'")
    return frozenset(ln.split()[0] for ln in lines[1:] if ln.strip())


def read_site_table(path: Path) -> list[SiteRecord]:
    records: list[SiteRecord] = []
    for i, ln in enumerate(Path(path).read_text().splitlines(), start=1):
        if not ln.strip():
            continue
        fields = ln.split("\t")
        if len(fields) != 3:
            raise SiteTableError(f"site table line {i}: expected 3 tab-separated fields")
        rsid, pos_s, af_s = fields
        try:
            pos = int(pos_s)
            af = float(af_s)
        except ValueError as e:
            raise SiteTableError(f"site table line {i}: {e}") from e
        if not (0.0 <= af <= 1.0):
            raise SiteTableError(f"site table line {i}: AF {af} outside [0, 1]")
        records.append(SiteRecord(rsid, pos, af))
    return records


def _apply_id_policy(records: list[SiteRecord]) -> tuple[list[SiteRecord], int, int]:
    """Policy steps 1-2: drop missing IDs, dedupe keep-first -> (usable, n_missing, n_dup)."""
    n_missing = 0
    n_dup = 0
    seen: set[str] = set()
    usable: list[SiteRecord] = []
    for r in records:
        if r.rsid in (".", ""):
            n_missing += 1
            continue
        if r.rsid in seen:
            n_dup += 1
            continue
        seen.add(r.rsid)
        usable.append(r)
    return usable, n_missing, n_dup


def partition_sites(
    records: list[SiteRecord],
    hm3_rsids: frozenset[str],
    maf_floor: float,
) -> SitePartition:
    usable, n_missing, n_dup = _apply_id_policy(records)

    array = tuple(r for r in usable if r.rsid in hm3_rsids)
    array_ids = {r.rsid for r in array}

    evaluation: list[SiteRecord] = []
    n_maf_excluded = 0
    for r in usable:
        if r.rsid in array_ids:
            continue
        if maf_floor <= r.af <= 1.0 - maf_floor:
            evaluation.append(r)
        else:
            n_maf_excluded += 1

    return SitePartition(
        array_sites=array,
        evaluation_sites=tuple(evaluation),
        n_missing_id=n_missing,
        n_duplicate_id=n_dup,
        n_eval_excluded_maf=n_maf_excluded,
    )


def partition_sites_by_thinning(
    records: list[SiteRecord],
    n_array: int,
    maf_floor: float,
) -> SitePartition:
    """Thin mode: n_array evenly spaced common sites become the array; remaining common
    sites are the evaluation set; rare sites (outside the MAF window) join neither
    (counted in n_eval_excluded_maf, which in this mode means 'excluded from both pools')."""
    usable, n_missing, n_dup = _apply_id_policy(records)

    common = sorted(
        (r for r in usable if maf_floor <= r.af <= 1.0 - maf_floor),
        key=lambda r: r.pos,
    )
    n_rare = len(usable) - len(common)

    if n_array < 2:
        raise SiteTableError(f"n_array must be >= 2, got {n_array}")
    if len(common) < 2 * n_array:
        raise SiteTableError(
            f"thinning needs >= 2x n_array common sites (evaluation must remain non-trivial); "
            f"have {len(common)} common sites for n_array={n_array}"
        )

    step = (len(common) - 1) / (n_array - 1)
    idx = sorted({round(i * step) for i in range(n_array)})
    if len(idx) != n_array:
        raise SiteTableError(
            f"even-spacing index collision: {len(idx)} unique of {n_array} requested"
        )
    idx_set = set(idx)
    array = tuple(common[i] for i in idx)
    evaluation = tuple(r for i, r in enumerate(common) if i not in idx_set)

    return SitePartition(
        array_sites=array,
        evaluation_sites=evaluation,
        n_missing_id=n_missing,
        n_duplicate_id=n_dup,
        n_eval_excluded_maf=n_rare,
    )


def write_partition(partition: SitePartition, outdir: Path) -> dict:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "array_sites.rsids.txt").write_text(
        "\n".join(r.rsid for r in partition.array_sites) + "\n"
    )
    (outdir / "evaluation_sites.rsids.txt").write_text(
        "\n".join(r.rsid for r in partition.evaluation_sites) + "\n"
    )
    summary = {
        "n_array_sites": len(partition.array_sites),
        "n_evaluation_sites": len(partition.evaluation_sites),
        "n_missing_id_excluded": partition.n_missing_id,
        "n_duplicate_id_excluded": partition.n_duplicate_id,
        "n_eval_excluded_by_maf": partition.n_eval_excluded_maf,
    }
    (outdir / "sites_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
