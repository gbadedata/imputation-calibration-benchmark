# imputation-calibration-benchmark

**Is Beagle's DR2 a calibrated confidence score - and does its calibration hold across
ancestries?** A Nextflow-wrapped benchmark on 1000 Genomes chromosome 20: simulate array
genotyping by masking WGS calls to a 20,000-site array-density site set, impute the masked targets against (A) an
ancestry-matched and (B) a deliberately mismatched (EUR panel, AFR targets) reference panel,
then grade every imputed genotype against the withheld WGS truth. DR2-vs-empirical-accuracy
calibration curves with bootstrap CIs, an abstention threshold chosen on Condition A and
transferred to Condition B, results reported as found.

Status: **Phase 1 (data derivation) - in progress.** Phase 0 verification is complete and
logged in [docs/verification.md](docs/verification.md).

## Known limitations (read first)

- **Within-cohort circularity (Condition A):** the matched reference panel and the targets
  are drawn from the same cohort (1000 Genomes), sharing sequencing platform, joint
  variant-calling pipeline, and cohort-level LD structure even after relatedness filtering.
  Condition A is therefore a favourable-case control, not an estimate of real-world
  matched-panel performance. The A-vs-B contrast is the claim; A's absolute numbers are not.
- **Single chromosome:** chr20 is a model system chosen for a 32 GB laptop budget.
  Calibration of per-variant confidence is measurable on one chromosome; nothing here is a
  polygenic-score result.
- **Array simulation:** the array is a deterministic uniform thinning of common
  (MAF >= 0.01) biallelic SNPs to 20,000 sites - array density and common-variant bias
  are mimicked; the site-selection strategy of real chips is not. HapMap3 intersection
  was the original design but is infeasible here: the 1kGP 20130502 v5b chr20 release
  carries no rsIDs (ID='.' on every row; docs/verification.md). The hm3 code path is
  retained and tested but unused.
- **Truth definition:** "truth" is the 1kGP phase 3 release call set, which has its own error
  rate; concordance with it is the measured quantity.

## Layout

    config/    Pydantic settings - every seed and threshold, ICB_* env-overridable
    src/       sample_design, sites, cli (Phase 1); accuracy/calibration/abstention arrive Phase 4
    scripts/   phase1.sh - the exact bcftools + CLI derivation pipeline
    tests/     synthetic-fixture tests, no network, no real data required
    docs/      verification.md - Phase 0 log
    data/      raw/ and derived/ are gitignored; PROVENANCE.md is the committed record

## Reproduce Phase 1

    python3 -m venv .venv && source .venv/bin/activate
    pip install -e ".[dev]"
    ruff check . && pytest
    # place raw inputs per data/PROVENANCE.md, then:
    bash scripts/phase1.sh
