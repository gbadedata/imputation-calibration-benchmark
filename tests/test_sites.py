"""Phase 1 tests: site partition policy on synthetic fixtures. No network, no real data."""

import pytest

from src.sites import (
    SiteRecord,
    SiteTableError,
    partition_sites,
    read_hm3_rsids,
    read_site_table,
    write_partition,
)

HM3 = frozenset({"rs1", "rs2", "rs3"})


def rec(rsid, pos, af):
    return SiteRecord(rsid, pos, af)


def test_read_hm3_parses_and_rejects_bad_header(tmp_path):
    good = tmp_path / "hm3.snplist"
    good.write_text("SNP A1 A2\nrs1 A G\nrs2 C T\n")
    assert read_hm3_rsids(good) == frozenset({"rs1", "rs2"})
    bad = tmp_path / "bad.snplist"
    bad.write_text("ID REF ALT\nrs1 A G\n")
    with pytest.raises(SiteTableError, match="header"):
        read_hm3_rsids(bad)


def test_read_site_table_parses_and_validates(tmp_path):
    path = tmp_path / "sites.tsv"
    path.write_text("rs1\t100\t0.25\nrs2\t200\t0.5\n")
    records = read_site_table(path)
    assert records == [rec("rs1", 100, 0.25), rec("rs2", 200, 0.5)]
    bad = tmp_path / "bad.tsv"
    bad.write_text("rs1\t100\t1.5\n")
    with pytest.raises(SiteTableError, match="outside"):
        read_site_table(bad)


def test_array_sites_are_hm3_members_only():
    records = [rec("rs1", 100, 0.3), rec("rs9", 200, 0.3)]
    p = partition_sites(records, HM3, maf_floor=0.01)
    assert [r.rsid for r in p.array_sites] == ["rs1"]
    assert [r.rsid for r in p.evaluation_sites] == ["rs9"]


def test_array_and_evaluation_disjoint():
    records = [rec(f"rs{i}", i * 100, 0.3) for i in range(1, 8)]
    p = partition_sites(records, HM3, maf_floor=0.01)
    array_ids = {r.rsid for r in p.array_sites}
    eval_ids = {r.rsid for r in p.evaluation_sites}
    assert not array_ids & eval_ids
    assert array_ids == {"rs1", "rs2", "rs3"}


def test_maf_floor_excludes_both_tails():
    records = [
        rec("rs10", 100, 0.005),  # below floor
        rec("rs11", 200, 0.01),  # exactly at floor - retained
        rec("rs12", 300, 0.99),  # exactly at ceiling - retained
        rec("rs13", 400, 0.995),  # above ceiling
        rec("rs14", 500, 0.5),
    ]
    p = partition_sites(records, HM3, maf_floor=0.01)
    kept = {r.rsid for r in p.evaluation_sites}
    assert kept == {"rs11", "rs12", "rs14"}
    assert p.n_eval_excluded_maf == 2


def test_maf_floor_not_applied_to_array_sites():
    # array sites keep rare variants: real chips genotype rare sites too,
    # and QC thresholds are Phase 2's job, applied via plink2 with logging
    records = [rec("rs1", 100, 0.001)]
    p = partition_sites(records, HM3, maf_floor=0.01)
    assert [r.rsid for r in p.array_sites] == ["rs1"]
    assert p.n_eval_excluded_maf == 0


def test_missing_ids_excluded_from_both_and_counted():
    records = [rec(".", 100, 0.3), rec("", 150, 0.3), rec("rs1", 200, 0.3)]
    p = partition_sites(records, HM3, maf_floor=0.01)
    assert p.n_missing_id == 2
    assert {r.rsid for r in p.array_sites} == {"rs1"}
    assert not p.evaluation_sites


def test_duplicate_rsids_keep_first_and_counted():
    records = [rec("rs9", 100, 0.3), rec("rs9", 999, 0.4), rec("rs1", 200, 0.3)]
    p = partition_sites(records, HM3, maf_floor=0.01)
    assert p.n_duplicate_id == 1
    eval_rs9 = [r for r in p.evaluation_sites if r.rsid == "rs9"]
    assert eval_rs9 == [rec("rs9", 100, 0.3)]  # first occurrence kept


def test_write_partition_summary_counts(tmp_path):
    records = [rec("rs1", 100, 0.3), rec("rs9", 200, 0.5), rec(".", 300, 0.5)]
    p = partition_sites(records, HM3, maf_floor=0.01)
    summary = write_partition(p, tmp_path)
    assert summary == {
        "n_array_sites": 1,
        "n_evaluation_sites": 1,
        "n_missing_id_excluded": 1,
        "n_duplicate_id_excluded": 0,
        "n_eval_excluded_by_maf": 0,
    }
    assert (tmp_path / "array_sites.rsids.txt").read_text() == "rs1\n"
    assert (tmp_path / "evaluation_sites.rsids.txt").read_text() == "rs9\n"


# ---- thin mode (primary for the v5b data, which carries no rsIDs) ----

from src.sites import partition_sites_by_thinning  # noqa: E402


def common_records(n, af=0.3, start_pos=1000):
    return [rec(f"20:{start_pos + i * 10}", start_pos + i * 10, af) for i in range(n)]


def test_thinning_selects_exact_count_and_disjoint():
    records = common_records(50)
    p = partition_sites_by_thinning(records, n_array=10, maf_floor=0.01)
    assert len(p.array_sites) == 10
    assert len(p.evaluation_sites) == 40
    array_ids = {r.rsid for r in p.array_sites}
    eval_ids = {r.rsid for r in p.evaluation_sites}
    assert not array_ids & eval_ids


def test_thinning_is_deterministic_and_spans_the_region():
    records = common_records(50)
    p1 = partition_sites_by_thinning(records, n_array=10, maf_floor=0.01)
    p2 = partition_sites_by_thinning(records, n_array=10, maf_floor=0.01)
    assert p1.array_sites == p2.array_sites
    positions = [r.pos for r in p1.array_sites]
    assert positions[0] == records[0].pos  # first common site included
    assert positions[-1] == records[-1].pos  # last common site included
    assert positions == sorted(positions)


def test_thinning_pool_is_common_only_rare_excluded_from_both():
    records = common_records(40) + [rec(f"20:{i}", i, 0.001) for i in range(1, 11)]
    p = partition_sites_by_thinning(records, n_array=10, maf_floor=0.01)
    assert p.n_eval_excluded_maf == 10  # in thin mode: rare sites excluded from both pools
    all_afs = [r.af for r in p.array_sites + p.evaluation_sites]
    assert all(0.01 <= af <= 0.99 for af in all_afs)


def test_thinning_applies_shared_id_policy():
    records = [rec(".", 5, 0.3), rec("20:100", 100, 0.3)] + common_records(40, start_pos=200)
    dup = rec("20:200", 999, 0.4)  # duplicate of first common_records id
    p = partition_sites_by_thinning(records + [dup], n_array=10, maf_floor=0.01)
    assert p.n_missing_id == 1
    assert p.n_duplicate_id == 1


def test_thinning_requires_headroom_over_n_array():
    with pytest.raises(SiteTableError, match="2x n_array"):
        partition_sites_by_thinning(common_records(15), n_array=10, maf_floor=0.01)
