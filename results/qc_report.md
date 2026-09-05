# Phase 2 QC report

Kinship: 33 samples removed at KING cutoff 0.0884 (2nd degree), across all AFR+EUR candidates, before the redraw.

Chip-cohort variant QC (targets only - see BUILD-SPEC 2.2 for why panels are not used in chip QC):

| Step | Threshold | Variants removed |
|---|---|---|
| Input | - | 20000 |
| Call rate (--geno) | 0.02 | 0 |
| MAF (--maf) | 0.01 | 2280 |
| HWE (--hwe) | 1e-06 | 22 |
| **Kept** | - | **17698** |

Samples in QC run: 100.

Accounting identity (input - removals == kept) verified at report build time.
