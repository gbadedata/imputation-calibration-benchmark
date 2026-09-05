"""BGEN round-trip verification (Phase 2.4).

Input: two files of `bcftools query -f '%POS[\\t%GT]\\n'` output - one from the original
VCF, one from the VCF re-exported after a BGEN round trip. Comparison is by per-site ALT
allele count, which is invariant to phase symbols ('|' vs '/') and heterozygote allele
order - the representation details a lossless round trip is allowed to change.

The verifier knows the failure mode Phase 0 caught on 2026-09-04: importing with the wrong
REF/ALT mode complements every count (c -> 2N - c, where N = sample count). It diagnoses
that pattern explicitly so the fix (ref-last, not ref-first) is named in the error.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


class RoundTripError(ValueError):
    """Raised when the round-tripped genotype content does not match the original."""


@dataclass(frozen=True)
class SiteCounts:
    pos: int
    alt_count: int
    n_samples: int


def alt_count_line(line: str) -> SiteCounts:
    """Parse one query line: POS<TAB>GT<TAB>GT... Missing alleles ('.') contribute 0."""
    fields = line.rstrip("\n").split("\t")
    pos = int(fields[0])
    count = 0
    for gt in fields[1:]:
        for allele in gt.replace("|", "/").split("/"):
            if allele == "1":
                count += 1
    return SiteCounts(pos=pos, alt_count=count, n_samples=len(fields) - 1)


def read_signature(path: Path) -> list[SiteCounts]:
    return [alt_count_line(ln) for ln in Path(path).read_text().splitlines() if ln.strip()]


def verify_round_trip(original: list[SiteCounts], round_tripped: list[SiteCounts]) -> str:
    if len(original) != len(round_tripped):
        raise RoundTripError(
            f"site count mismatch: {len(original)} original vs {len(round_tripped)} "
            "round-tripped - variants were dropped or duplicated"
        )
    if not original:
        raise RoundTripError("empty signature files - nothing was compared")

    mismatches = [
        (a, b) for a, b in zip(original, round_tripped, strict=True) if a != b
    ]
    if not mismatches:
        return f"ROUND TRIP CLEAN: {len(original)} sites, genotype content identical"

    # Diagnose the known allele-order failure: every count complemented to 2N - c
    if all(
        a.pos == b.pos and b.alt_count == 2 * a.n_samples - a.alt_count
        for a, b in zip(original, round_tripped, strict=True)
    ):
        raise RoundTripError(
            "ALLELE-ORDER FLIP: every alt count is complemented (c -> 2N - c). "
            "The BGEN import REF/ALT mode is wrong - plink2 a7.4 bgen-1.2 export is "
            "ref-LAST; import must use 'ref-last' (docs/verification.md V3, 2026-09-04)."
        )
    a, b = mismatches[0]
    raise RoundTripError(
        f"{len(mismatches)} of {len(original)} sites differ; first at POS {a.pos}: "
        f"original alt={a.alt_count}, round-tripped alt={b.alt_count}. Not a uniform "
        "complement - investigate dosage handling."
    )
