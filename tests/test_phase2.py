"""Phase 2 tests: QC report machinery and BGEN round-trip verifier. Synthetic inputs only."""

import pytest

from src.bgen_roundtrip import (
    RoundTripError,
    SiteCounts,
    alt_count_line,
    verify_round_trip,
)
from src.qc_report import (
    QCReportError,
    build_report,
    count_id_file,
    count_snplist,
    parse_plink_log,
    render_markdown,
)

SYNTH_LOG = """PLINK v2.0.0-a.7.4LM AVX2 Intel (18 Aug 2026)
--vcf: 20000 variants scanned.
100 samples (0 females, 0 males, 100 ambiguous; 100 founders) loaded from
20000 variants loaded from data/derived/qc/plink2-temporary.pvar.
--geno: 3 variants removed due to missing genotype data.
1200 variants removed due to allele frequency threshold(s)
--hwe: 7 variants removed due to Hardy-Weinberg exact test (founders only).
18790 variants remaining after main filters.
"""


def test_parse_plink_log_extracts_all_counts():
    c = parse_plink_log(SYNTH_LOG)
    assert c == {
        "n_input_variants": 20000,
        "n_samples": 100,
        "n_geno_removed": 3,
        "n_maf_removed": 1200,
        "n_hwe_removed": 7,
    }


def test_parse_plink_log_missing_removal_lines_default_zero():
    log = "--vcf: 500 variants scanned.\n60 samples (0 females) loaded\n"
    c = parse_plink_log(log)
    assert c["n_geno_removed"] == 0
    assert c["n_maf_removed"] == 0
    assert c["n_hwe_removed"] == 0


def test_parse_plink_log_missing_input_count_raises():
    with pytest.raises(QCReportError, match="n_input_variants"):
        parse_plink_log("no useful content here\n")


def test_count_id_and_snplist_files(tmp_path):
    idf = tmp_path / "king.cutoff.out.id"
    idf.write_text("#FID\tIID\n0\tNA001\n0\tNA002\n")
    assert count_id_file(idf) == 2
    snp = tmp_path / "kept.snplist"
    snp.write_text("20:100\n20:200\n20:300\n")
    assert count_snplist(snp) == 3


def make_report_files(tmp_path, kept_lines):
    snp = tmp_path / "kept.snplist"
    snp.write_text("".join(f"20:{i}\n" for i in range(kept_lines)))
    idf = tmp_path / "king.id"
    idf.write_text("#IID\nNA001\n")
    return snp, idf


def test_build_report_accounting_identity_holds(tmp_path):
    snp, idf = make_report_files(tmp_path, kept_lines=18790)
    r = build_report(
        SYNTH_LOG, snplist_path=snp, king_id_path=idf,
        geno=0.02, maf=0.01, hwe=1e-6, king_cutoff=0.0884,
    )
    assert r.n_kept_variants == 18790
    assert r.n_kinship_removed_samples == 1


def test_build_report_accounting_identity_violation_raises(tmp_path):
    snp, idf = make_report_files(tmp_path, kept_lines=18000)  # inconsistent with log
    with pytest.raises((QCReportError, Exception), match="accounting"):
        build_report(
            SYNTH_LOG, snplist_path=snp, king_id_path=idf,
            geno=0.02, maf=0.01, hwe=1e-6, king_cutoff=0.0884,
        )


def test_render_markdown_contains_all_numbers(tmp_path):
    snp, idf = make_report_files(tmp_path, kept_lines=18790)
    r = build_report(
        SYNTH_LOG, snplist_path=snp, king_id_path=idf,
        geno=0.02, maf=0.01, hwe=1e-6, king_cutoff=0.0884,
    )
    md = render_markdown(r)
    for token in ("20000", "1200", "18790", "0.0884"):
        assert token in md


# ---- BGEN round-trip verifier ----


def test_alt_count_line_phase_and_order_invariant():
    assert alt_count_line("100\t0|1\t1|0\t1/1").alt_count == 4
    assert alt_count_line("100\t1|0\t0/1\t1|1").alt_count == 4


def test_alt_count_line_missing_alleles_contribute_zero():
    sc = alt_count_line("100\t./.\t0|1\t.")
    assert sc.alt_count == 1
    assert sc.n_samples == 3


def test_verify_round_trip_clean():
    sig = [SiteCounts(100, 4, 3), SiteCounts(200, 0, 3)]
    verdict = verify_round_trip(sig, list(sig))
    assert "CLEAN" in verdict and "2 sites" in verdict


def test_verify_round_trip_diagnoses_allele_order_flip():
    original = [SiteCounts(100, 1, 3), SiteCounts(200, 5, 3)]
    flipped = [SiteCounts(100, 5, 3), SiteCounts(200, 1, 3)]  # c -> 2N - c with N=3
    with pytest.raises(RoundTripError, match="ref-last"):
        verify_round_trip(original, flipped)


def test_verify_round_trip_reports_non_complement_mismatch():
    original = [SiteCounts(100, 1, 3), SiteCounts(200, 5, 3)]
    corrupt = [SiteCounts(100, 2, 3), SiteCounts(200, 5, 3)]
    with pytest.raises(RoundTripError, match="POS 100"):
        verify_round_trip(original, corrupt)


def test_verify_round_trip_site_count_mismatch():
    with pytest.raises(RoundTripError, match="site count"):
        verify_round_trip([SiteCounts(100, 1, 3)], [])
