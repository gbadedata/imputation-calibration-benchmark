"""Phase 1 tests: sample design invariants on synthetic panel fixtures. No network, no real data."""

import pytest

from config.settings import Settings
from src.sample_design import (
    SampleDesignError,
    design_samples,
    read_panel,
    read_sample_list,
    write_design,
)


def small_settings(**overrides) -> Settings:
    base = {"n_targets": 3, "seed": 42}
    base.update(overrides)
    return Settings(**base)


def make_panel_file(tmp_path, rows):
    lines = ["sample\tpop\tsuper_pop\tgender"]
    lines += [f"{s}\t{p}\t{sp}\tmale" for s, p, sp in rows]
    path = tmp_path / "panel.txt"
    path.write_text("\n".join(lines) + "\n")
    return path


AFR_ROWS = [(f"AFR{i:03d}", "YRI", "AFR") for i in range(8)]
EUR_ROWS = [(f"EUR{i:03d}", "GBR", "EUR") for i in range(5)]
EAS_ROWS = [(f"EAS{i:03d}", "CHB", "EAS") for i in range(2)]


def test_read_panel_parses_rows_and_superpops(tmp_path):
    path = make_panel_file(tmp_path, AFR_ROWS + EUR_ROWS + EAS_ROWS)
    panel = read_panel(path)
    assert len(panel) == 15
    assert sum(1 for s in panel if s.super_pop == "AFR") == 8
    assert sum(1 for s in panel if s.super_pop == "EUR") == 5


def test_read_panel_rejects_bad_header(tmp_path):
    path = tmp_path / "panel.txt"
    path.write_text("id\tgroup\tregion\nA\tYRI\tAFR\n")
    with pytest.raises(SampleDesignError, match="header"):
        read_panel(path)


def test_read_panel_rejects_duplicate_sample(tmp_path):
    path = make_panel_file(tmp_path, [("S1", "YRI", "AFR"), ("S1", "YRI", "AFR")])
    with pytest.raises(SampleDesignError, match="duplicate"):
        read_panel(path)


def test_target_draw_reproducible_from_seed(tmp_path):
    panel = read_panel(make_panel_file(tmp_path, AFR_ROWS + EUR_ROWS))
    d1 = design_samples(panel, set(), small_settings())
    d2 = design_samples(panel, set(), small_settings())
    assert d1.targets == d2.targets
    d3 = design_samples(panel, set(), small_settings(seed=7))
    assert d3.targets != d1.targets  # 8-choose-3 = 56 draws; seed collision effectively impossible


def test_draw_independent_of_panel_row_order(tmp_path):
    panel_fwd = read_panel(make_panel_file(tmp_path, AFR_ROWS + EUR_ROWS))
    shuffled = list(reversed(AFR_ROWS)) + EUR_ROWS
    (tmp_path / "panel.txt").unlink()
    panel_rev = read_panel(make_panel_file(tmp_path, shuffled))
    s = small_settings()
    t_fwd = design_samples(panel_fwd, set(), s).targets
    t_rev = design_samples(panel_rev, set(), s).targets
    assert t_fwd == t_rev


def test_counts_superpop_purity_and_disjointness(tmp_path):
    panel = read_panel(make_panel_file(tmp_path, AFR_ROWS + EUR_ROWS + EAS_ROWS))
    d = design_samples(panel, set(), small_settings())
    assert len(d.targets) == 3
    assert len(d.panel_a) == 5  # 8 AFR - 3 targets
    assert len(d.panel_b) == 5  # all EUR
    assert all(t.startswith("AFR") for t in d.targets)
    assert all(a.startswith("AFR") for a in d.panel_a)
    assert all(b.startswith("EUR") for b in d.panel_b)
    assert not set(d.targets) & set(d.panel_a)
    assert not set(d.targets) & set(d.panel_b)
    assert not set(d.panel_a) & set(d.panel_b)


def test_kinship_removals_excluded_everywhere(tmp_path):
    panel = read_panel(make_panel_file(tmp_path, AFR_ROWS + EUR_ROWS))
    removed = {"AFR000", "EUR000"}
    d = design_samples(panel, removed, small_settings())
    all_assigned = set(d.targets) | set(d.panel_a) | set(d.panel_b)
    assert not all_assigned & removed
    assert d.n_removed_by_kinship == 2


def test_insufficient_target_pool_raises(tmp_path):
    panel = read_panel(make_panel_file(tmp_path, AFR_ROWS[:3] + EUR_ROWS))
    with pytest.raises(SampleDesignError, match="eligible AFR"):
        design_samples(panel, set(), small_settings())  # needs > 3 AFR, has exactly 3


def test_empty_panel_b_raises(tmp_path):
    panel = read_panel(make_panel_file(tmp_path, AFR_ROWS))
    with pytest.raises(SampleDesignError, match="panel B"):
        design_samples(panel, set(), small_settings())


def test_write_design_round_trips(tmp_path):
    panel = read_panel(make_panel_file(tmp_path, AFR_ROWS + EUR_ROWS))
    s = small_settings()
    d = design_samples(panel, set(), s)
    summary = write_design(d, tmp_path / "out", s)
    assert summary["n_targets"] == 3
    assert summary["seed"] == 42
    targets_back = read_sample_list(tmp_path / "out" / "targets_afr.txt")
    assert targets_back == set(d.targets)
    panel_a_back = read_sample_list(tmp_path / "out" / "panel_a_afr.txt")
    assert panel_a_back == set(d.panel_a)
