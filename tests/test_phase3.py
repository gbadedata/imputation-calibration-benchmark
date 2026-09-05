"""Phase 3 tests: imputation output validation on synthetic fixtures."""

import pytest

from config.settings import Settings
from src.imputation_check import (
    ImputationCheckError,
    check_condition,
    header_declares,
    parse_dr2_table,
)

HEADER_OK = (
    '##INFO=<ID=DR2,Number=A,Type=Float,Description="Dosage R-Squared">\n'
    '##FORMAT=<ID=DS,Number=A,Type=Float,Description="estimated ALT dose">\n'
)
TARGETS = ["T1", "T2", "T3"]
EVAL = {"20:100", "20:200", "20:300", "20:400"}


def base_kwargs(**over):
    kw = dict(
        condition="matched",
        header_text=HEADER_OK,
        sample_ids=list(TARGETS),
        target_ids=list(TARGETS),
        output_marker_ids=EVAL | {"20:999"},
        eval_site_ids=set(EVAL),
        dr2_values={"20:100": 0.9, "20:200": 0.5, "20:300": 0.99, "20:400": 0.1},
        n_dr2_missing=0,
        settings=Settings(),
    )
    kw.update(over)
    return kw


def test_header_declares_detects_fields():
    assert header_declares(HEADER_OK, "DR2")
    assert header_declares(HEADER_OK, "DS")
    assert not header_declares(HEADER_OK, "GP")


def test_parse_dr2_table_values_and_missing(tmp_path):
    p = tmp_path / "dr2.tsv"
    p.write_text("20:100\t0.97\n20:200\t.\n20:300\t0.12\n")
    values, n_missing = parse_dr2_table(p)
    assert values == {"20:100": 0.97, "20:300": 0.12}
    assert n_missing == 1


def test_parse_dr2_table_rejects_multiallelic(tmp_path):
    p = tmp_path / "dr2.tsv"
    p.write_text("20:100\t0.9,0.8\n")
    with pytest.raises(ImputationCheckError, match="biallelic"):
        parse_dr2_table(p)


def test_matched_condition_passes_clean_fixture():
    m = check_condition(**base_kwargs())
    assert m.eval_coverage == 1.0
    assert m.samples_match_targets
    assert m.n_dr2_out_of_range == 0


def test_missing_dr2_declaration_fails_both_conditions():
    bad_header = '##FORMAT=<ID=DS,Number=A,Type=Float,Description="x">\n'
    for cond in ("matched", "mismatched"):
        with pytest.raises(ImputationCheckError, match="DR2 not declared"):
            check_condition(**base_kwargs(condition=cond, header_text=bad_header))


def test_sample_mismatch_fails():
    with pytest.raises(ImputationCheckError, match="samples"):
        check_condition(**base_kwargs(sample_ids=["T1", "T2", "OTHER"]))


def test_dr2_out_of_range_fails():
    with pytest.raises(ImputationCheckError, match="outside"):
        check_condition(**base_kwargs(dr2_values={"20:100": 1.2}))


def test_matched_coverage_floor_enforced():
    # only 2 of 4 eval sites present -> 50% < 95% floor
    with pytest.raises(ImputationCheckError, match="matched coverage"):
        check_condition(**base_kwargs(output_marker_ids={"20:100", "20:200"}))


def test_mismatched_condition_reports_low_coverage_without_failing():
    # 75% coverage: fails matched, passes mismatched (floor 0.5) and is reported
    kw = base_kwargs(
        condition="mismatched",
        output_marker_ids={"20:100", "20:200", "20:300"},
        dr2_values={"20:100": 0.9, "20:200": 0.5, "20:300": 0.99},
    )
    m = check_condition(**kw)
    assert m.eval_coverage == 0.75
    with pytest.raises(ImputationCheckError):
        check_condition(**base_kwargs(output_marker_ids={"20:100", "20:200", "20:300"}))


def test_mismatched_catastrophic_floor_still_fails():
    with pytest.raises(ImputationCheckError, match="catastrophic"):
        check_condition(
            **base_kwargs(condition="mismatched", output_marker_ids={"20:100"})
        )


def test_matched_dr2_missing_fraction_gated():
    with pytest.raises(ImputationCheckError, match="DR2-missing"):
        check_condition(
            **base_kwargs(dr2_values={"20:100": 0.9, "20:200": 0.5}, n_dr2_missing=2)
        )
