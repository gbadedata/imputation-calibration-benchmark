"""Settings behave as the single source of truth and honour the ICB_ env prefix."""

from config.settings import Settings


def test_defaults_match_build_spec():
    s = Settings()
    assert s.seed == 42
    assert s.n_targets == 100
    assert s.target_superpop == "AFR"
    assert s.mismatched_panel_superpop == "EUR"
    assert s.eval_maf_floor == 0.01
    assert s.king_cutoff == 0.0884
    assert s.dr2_bins == 10
    assert s.bootstrap_reps == 1000


def test_env_prefix_overrides(monkeypatch):
    monkeypatch.setenv("ICB_SEED", "7")
    monkeypatch.setenv("ICB_N_TARGETS", "10")
    s = Settings()
    assert s.seed == 7
    assert s.n_targets == 10
