"""Site design: partition chr20 biallelic SNPs into array sites and evaluation sites.

Inputs (both produced by scripts/phase1.sh, never parsed from VCF here):
- w_hm3.snplist: header 'SNP A1 A2', whitespace-separated, rsIDs genome-wide
- site table TSV: ID<TAB>POS<TAB>AF per biallelic chr20 SNP (bcftools query output)

Policy, in order (each exclusion is counted, never silent - BUILD-SPEC Phase 1.2):
1. rows with missing IDs ('.' or empty) are excluded from BOTH site sets: without an
   rsID a site cannot be matched to the hm3 list nor reliably re-joined downstream
2. duplicate rsIDs keep the first occurrence; later occurrences are excluded from both
3. array sites   = remaining sites whose rsID is in the hm3 list
4. evaluation    = remaining sites NOT in the array set with maf_floor <= AF <= 1 - maf_floor

Invariant: array_sites and evaluation_sites are disjoint by construction; tests assert it.
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


def partition_sites(
    records: list[SiteRecord],
    hm3_rsids: frozenset[str],
    maf_floor: float,
) -> SitePartition:
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
