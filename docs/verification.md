# Phase 0 verification log

Date: 2026-09-04. Machine: OJ-ODIMAYO, WSL2 Ubuntu 24.04, 32 GB RAM, i7-1165G7.

## V1 - 1000 Genomes chr20 + panel: PASS

- Directory: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
- chr20 genotypes file is **v5b**: `ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`.
  Community scripts (MRCIEU get_1000g.sh and others, ~2017) cite v5a for chr20; the live
  server 404s v5a and serves v5b. Resolved by listing the directory itself
  (`wget -qO- <dir>/ | grep -o 'ALL\.chr20[^"<]*'`). Lesson: filenames are verified against
  the live server, not tutorials.
- Downloaded .tbi predates index count metadata; working index regenerated locally with
  bcftools (`bcftools index -f -t`). Downloaded .tbi checksum retained in CHECKSUMS.txt.
- `bcftools index -s`: chr20, length 63025520, **1,812,841 variants** (expected 1.7-1.9M).
- Panel `integrated_call_samples_v3.20130502.ALL.panel`: md5 7ee5675553088230530a7fe88c22f201,
  2,504 samples; super_pop counts AFR 661, AMR 347, EAS 504, EUR 503, SAS 489.

## V2 - HapMap3 site list: PASS

- Source: Zenodo record 7773502 (KCL archival copy of Alkes Group w_hm3.snplist; DOI-backed).
- md5 of w_hm3.snplist.gz: 153ecc2bcfa740afafe656e6a384d769 (matches Zenodo's published md5).
- 1,217,311 SNPs; header `SNP A1 A2`; rsIDs only, no positions -> array-site intersection is
  by VCF ID column; genome-wide list self-selects chr20 against a chr20-only VCF.

## V3 - plink2 + BGEN round trip: PENDING PASTE-BACK

- Build: alpha 7.4 (18 Aug 2026), Linux AVX2 Intel,
  https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_avx2_20260818.zip
- Flags --geno/--maf/--hwe/--king-cutoff/--export bgen-1.2 confirmed in current docs.
- BGEN import requires explicit REF/ALT mode (`ref-first`).
- TODO: paste plink2 --version string and the round-trip PASS line, then flip to PASS.

## V4 - Beagle 5.5 dry run: PASS

- Jar: beagle.27Feb25.75f.jar (version 5.5), 316,844 bytes, from faculty.washington.edu.
- Genetic map: plink.GRCh37.map.zip from bochet.gcc.biostat.washington.edu (23,913,120 bytes);
  using plink.chr20.GRCh37.map.
- Java: OpenJDK 21.0.12 (requirement is 1.8+).
- Dry run (20 targets, 200-sample ref, every-5th-site thinning, -Xmx8g, nthreads=6):
  13,367 reference markers, 2,673 study markers, total time 1 s,
  peak RSS **271,152 kB** (~271 MB) at toy scale. Phase 3 has its own full-scale dry-run gate.
- Output header declares DR2 (INFO) and DS (FORMAT).
- **DR2 is `Number=A`** (per-ALT-allele, not per-site). Biallelic-only design collapses this
  to one value per variant, but Phase 4 parsing must index DR2 explicitly, with a test.
- Beagle 5.5 requires the ref panel fully phased with '|' separators and no missing genotypes;
  1kGP phase3 satisfies this, and all Phase 1 subsetting is via bcftools view (phase-preserving).
