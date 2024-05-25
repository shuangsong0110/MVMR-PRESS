# MVMR-PRESS
Implementation codes for "A robust penalized-regression-based method for multivariable Mendelian randomization using GWAS summary statistics"

## Table of contents
* [Prerequisites](#white_check_mark-prerequisites)
* [Installation](#hammer_and_wrench-installation)
* [Prepare GWAS summary statistics](#scroll-prepare-gwas-summary-statistics)

## :white_check_mark: Prerequisites

The software is developed and tested in Linux and Windows environments.

- R (>=3.6)
- GNU Scientific Library (GSL) (>=2.3)

## :hammer_and_wrench: Installation

```r
devtools::install_github("shuangsong0110/MVMR-PRESS")
```

## :scroll: Prepare GWAS summary statistics
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
