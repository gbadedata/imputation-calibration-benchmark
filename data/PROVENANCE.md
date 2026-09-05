# Data provenance

All raw inputs live in `data/raw/` (gitignored); this file is the committed record.
Derived-file section is appended by hand after each `scripts/phase1.sh` run.

## Raw inputs (downloaded 2026-09-04)

| File | Source | Checksum / size |
|---|---|---|
| ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ | see data/raw/CHECKSUMS.txt; 1,812,841 variants, 2,504 samples |
| (same).tbi | same directory | upstream index lacks count metadata; working index regenerated locally, bcftools `index -f -t` |
| integrated_call_samples_v3.20130502.ALL.panel | same directory | md5 7ee5675553088230530a7fe88c22f201; 2,504 samples |
| w_hm3.snplist.gz | Zenodo record 7773502 (KCL archival copy of Alkes Group list), https://zenodo.org/records/7773502 | md5 153ecc2bcfa740afafe656e6a384d769; 1,217,311 SNPs |
| beagle.27Feb25.75f.jar (v5.5) | https://faculty.washington.edu/browning/beagle/ | 316,844 bytes |
| plink.GRCh37.map.zip | http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/ | 23,913,120 bytes; using plink.chr20.GRCh37.map |
| plink2 alpha 7.4 (18 Aug 2026), Linux AVX2 | https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_avx2_20260818.zip | version string in docs/verification.md |

Filename note: chr20 genotypes are **v5b** on the live server; widely-circulated scripts cite
v5a and 404. Verified by directory listing on 2026-09-04 (docs/verification.md, V1).

Tool versions: bcftools (record `bcftools --version | head -1` here after V3 paste-back),
OpenJDK 21.0.12, Beagle 5.5 (27Feb25.75f), plink2 alpha 7.4.

## Derived files (append per run)

<!-- After each scripts/phase1.sh run, append: date, git commit of the code that ran,
     samples_summary.json and sites_summary.json contents, md5sums of the four derived
     VCFs, and whether a kinship-removal list was applied. -->

### Derivation run 2026-09-05, code commit 4c918e9
```
{
  "seed": 42,
  "n_targets": 100,
  "n_panel_a": 561,
  "n_panel_b": 503,
  "n_removed_by_kinship": 0,
  "target_superpop": "AFR",
  "panel_b_superpop": "EUR"
}
{
  "n_array_sites": 20000,
  "n_evaluation_sites": 246037,
  "n_missing_id_excluded": 0,
  "n_duplicate_id_excluded": 1,
  "n_eval_excluded_by_maf": 1473277
}
97bbd094630f414bf7eddf281a5429d7  data/derived/truth_targets.vcf.gz
3dd135d39a259ab2f41a3652c28c1cba  data/derived/array_targets.vcf.gz
2973186f5c798496e49b0cbb0691078c  data/derived/panelA.vcf.gz
38fc6b76a211fb49071dfca896c2b65e  data/derived/panelB.vcf.gz
```
