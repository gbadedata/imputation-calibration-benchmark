"""Sample design: draw imputation targets and reference panels from the 1kGP panel file.

Design (see BUILD-SPEC Phase 1.1):
- targets: n_targets samples of target_superpop (AFR), seeded random draw
- panel A (matched): all remaining eligible samples of matched_panel_superpop
- panel B (mismatched): all eligible samples of mismatched_panel_superpop
- 'eligible' = present in the panel file and not in the kinship-removal list
  (the removal list arrives from Phase 2.2; pass an empty set on the first run
  and re-run after --king-cutoff finalises removals)

Hard invariants, enforced here and again by tests:
- targets, panel A, panel B are pairwise disjoint
- panels are superpopulation-pure
- the draw is reproducible from the seed alone
"""

from __future__ import annotations

import json
import random
from dataclasses import dataclass
from pathlib import Path

from config.settings import Settings

PANEL_HEADER_PREFIX = ("sample", "pop", "super_pop")


@dataclass(frozen=True)
class PanelSample:
    sample: str
    pop: str
    super_pop: str


@dataclass(frozen=True)
class SampleDesign:
    targets: tuple[str, ...]
    panel_a: tuple[str, ...]
    panel_b: tuple[str, ...]
    n_removed_by_kinship: int


class SampleDesignError(ValueError):
    """Raised when the panel file or the resulting design violates an invariant."""


def read_panel(path: Path) -> list[PanelSample]:
    """Parse integrated_call_samples_v3-style panel file (TSV, header row)."""
    lines = Path(path).read_text().splitlines()
    if not lines:
        raise SampleDesignError(f"panel file is empty: {path}")
    header = tuple(lines[0].split("\t")[:3])
    if header != PANEL_HEADER_PREFIX:
        raise SampleDesignError(
            f"unexpected panel header {header!r}; expected columns {PANEL_HEADER_PREFIX!r}"
        )
    samples: list[PanelSample] = []
    seen: set[str] = set()
    for i, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            raise SampleDesignError(f"panel line {i} has {len(fields)} fields, expected >= 3")
        sid, pop, superpop = fields[0], fields[1], fields[2]
        if sid in seen:
            raise SampleDesignError(f"duplicate sample id in panel: {sid}")
        seen.add(sid)
        samples.append(PanelSample(sid, pop, superpop))
    return samples


def read_sample_list(path: Path) -> set[str]:
    """One sample id per line; blank lines ignored."""
    return {ln.strip() for ln in Path(path).read_text().splitlines() if ln.strip()}


def design_samples(
    panel: list[PanelSample],
    removals: set[str],
    settings: Settings,
) -> SampleDesign:
    eligible = [s for s in panel if s.sample not in removals]
    n_removed = len(panel) - len(eligible)

    # sorted() before sampling makes the draw independent of panel-file row order
    matched_pool = sorted(
        s.sample for s in eligible if s.super_pop == settings.target_superpop
    )
    if settings.matched_panel_superpop != settings.target_superpop:
        raise SampleDesignError(
            "design assumes panel A shares the target superpopulation; "
            f"got target={settings.target_superpop} panelA={settings.matched_panel_superpop}"
        )
    if len(matched_pool) < settings.n_targets + 1:
        raise SampleDesignError(
            f"need > {settings.n_targets} eligible {settings.target_superpop} samples "
            f"(targets + non-empty panel A); have {len(matched_pool)}"
        )

    rng = random.Random(settings.seed)
    targets = tuple(sorted(rng.sample(matched_pool, settings.n_targets)))
    panel_a = tuple(sorted(set(matched_pool) - set(targets)))
    panel_b = tuple(
        sorted(
            s.sample
            for s in eligible
            if s.super_pop == settings.mismatched_panel_superpop
        )
    )
    if not panel_b:
        raise SampleDesignError(
            f"no eligible samples for panel B ({settings.mismatched_panel_superpop})"
        )

    design = SampleDesign(targets, panel_a, panel_b, n_removed)
    _validate(design)
    return design


def _validate(d: SampleDesign) -> None:
    t, a, b = set(d.targets), set(d.panel_a), set(d.panel_b)
    if t & a:
        raise SampleDesignError(f"targets overlap panel A: {sorted(t & a)[:5]}")
    if t & b:
        raise SampleDesignError(f"targets overlap panel B: {sorted(t & b)[:5]}")
    if a & b:
        raise SampleDesignError(f"panel A overlaps panel B: {sorted(a & b)[:5]}")


def write_design(design: SampleDesign, outdir: Path, settings: Settings) -> dict:
    """Write sample lists + a machine-readable summary; return the summary dict."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "targets_afr.txt").write_text("\n".join(design.targets) + "\n")
    (outdir / "panel_a_afr.txt").write_text("\n".join(design.panel_a) + "\n")
    (outdir / "panel_b_eur.txt").write_text("\n".join(design.panel_b) + "\n")
    summary = {
        "seed": settings.seed,
        "n_targets": len(design.targets),
        "n_panel_a": len(design.panel_a),
        "n_panel_b": len(design.panel_b),
        "n_removed_by_kinship": design.n_removed_by_kinship,
        "target_superpop": settings.target_superpop,
        "panel_b_superpop": settings.mismatched_panel_superpop,
    }
    (outdir / "samples_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
