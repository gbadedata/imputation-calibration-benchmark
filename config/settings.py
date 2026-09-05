"""Central configuration. Every seed and threshold lives here, overridable via ICB_* env vars.

Example: ICB_SEED=7 python -m src.cli design-samples ...
"""

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_prefix="ICB_")

    # Sample design (Phase 1)
    seed: int = 42
    n_targets: int = 100
    target_superpop: str = "AFR"
    matched_panel_superpop: str = "AFR"  # Condition A reference
    mismatched_panel_superpop: str = "EUR"  # Condition B reference

    # Site design (Phase 1)
    eval_maf_floor: float = 0.01  # evaluation sites require floor <= AF <= 1 - floor
    n_array_sites: int = 20000  # thinning-mode array size (v5b release carries no rsIDs)

    # QC (Phase 2) - recorded here now so the report generator has one source of truth
    geno_missing_max: float = 0.02
    maf_min: float = 0.01
    hwe_p: float = 1e-6
    king_cutoff: float = 0.0884  # 2nd-degree threshold

    # Calibration (Phase 4)
    dr2_bins: int = 10
    bootstrap_reps: int = 1000
    abstention_accuracy_floor: float = 0.8  # CI lower bound required to retain a DR2 bin
