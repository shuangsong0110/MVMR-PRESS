# MVMR-PRESS
Implementation codes for "A robust penalized-regression-based method for multivariable Mendelian randomization using GWAS summary statistics"

## Table of contents
* [Prerequisites](#white_check_mark-prerequisites)
* [Installation](#hammer_and_wrench-installation)
* [Run MVMR-PRESS](#rocket-run-mvmr-press)
* [Output](#key-output)
* [GWAS summary statistics format](#scroll-gwas-summary-statistics-format)
* [A toy example](#bulb-a-toy-example)

## :white_check_mark: Prerequisites

The software is developed and tested in Linux and Windows environments.

- R (>=3.6)
- GNU Scientific Library (GSL) (>=2.3)

## :hammer_and_wrench: Installation

```r
devtools::install_github("shuangsong0110/MVMR-PRESS")
```

## :rocket: Run MVMR-PRESS
```
library(MVMRPRESS)
res <- run_mvmrpress(summs_exposure = SUMMARY_STATS_FOR_EXPOSURE (required),
                     summs_outcome = SUMMARY_STATS_FOR_OUTCOME (required),
                     snp.use = INSTRUMENTAL_VARIABLES (required),
                     P_hat = CORRELATION_MATRIX_OF_EXPOSURES (required),
                     para_theta=PARAMETER_THETA,
                     para_gamma=PARAMETER_GAMMA,
                     n_boots=BOOTSTRAP_NUMBER,
                     n_cores=NCORES)

```

- SUMMARY_STATS_FOR_EXPOSURE (required): GWAS summary statistics for exposures, should be a list, each element is a data.frame with headers SNP, A1, A2, BETA, SE (see [GWAS summary statistics format](#scroll-gwas-summary-statistics-format)).


- SUMMARY_STATS_FOR_OUTCOME (required): GWAS summary statistics for the outcome, should be a data.frame with headers SNP, A1, A2, BETA, SE (see [GWAS summary statistics format](#scroll-gwas-summary-statistics-format)).


- INSTRUMENTAL_VARIABLES (required): rsid of the IVs.

- CORRELATION_MATRIX_OF_EXPOSURES (required): correlation matrix of exposures which can be estimated with LDSC. Should be an identity matrix if no samples overlapping. 

- PARAMETER_THETA (optional): sparsity parameter for theta. Default=0.01.

- PARAMETER_GAMMA (optional): sparsity parameter for gamma. Default=10.
- BOOTSTRAP_NUMBER (optional): number of bootstrap. Default=100.
- NCORES (optional): number of cores run in parallel. Default=100.

## :key: Output
The `run_mvmrpress` returns a list with 4 elements:
- beta_est: estimated effect sizes for all exposures
- beta_est_SE: estimated SE of effect sizes
- beta_Z: Z scores for all exposures
- beta_P: P values for all exposures (adjusted for multiple testing)


## :scroll: GWAS summary statistics format
Please prepare the GWAS summary statistics for the outcome and each of the exposure  in the following format (including the header line):
```
   SNP        A1   A2       BETA         SE         
 rs4040617    G     A      -0.199       0.12
 rs4075116    C     T       0.646       0.05
 rs9442385    T     G      -0.016       0.37
    ...
```
**SNP**: SNP rsid

**A1**: reference allele

**A2**: alternative allele

**BETA**: GWAS effect sizes

**SE**: Standard errors of GWAS effect sizes

##  :bulb: A toy example
