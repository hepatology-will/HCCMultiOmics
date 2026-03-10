
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HCCMultiOmics

<!-- badges: start -->

[![R-CMD-check](https://github.com/hepatology-will/HCCMultiOmics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hepatology-will/HCCMultiOmics/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

HCCMultiOmics is an integrative AutoML framework for multiscale
biological discovery in Hepatocellular Carcinoma (HCC). The package
provides a comprehensive and automated workflow designed to translate
multi-omics data into actionable biological mechanisms.

## Features

- **Module A: Bulk Analysis**: Ensemble machine learning engine (LASSO,
  RSF, XGBoost) for robust prognostic biomarker discovery
- **Module B: Single-cell Analysis**: Adaptive single-cell mechanistic
  decoding module to profile intracellular trajectories, stemness, and
  intercellular communication
- **Integrated Workflow**: Seamless transition from bulk to single-cell
  analysis
- **Automated Visualization**: Comprehensive plotting functions for
  results interpretation

## Installation

You can install the development version of HCCMultiOmics from
[GitHub](https://github.com/hepatology-will/HCCMultiOmics) with:

``` r
# install.packages("devtools")
devtools::install_github("yourusername/HCCMultiOmics")
```

## Data Requirements

For SCENIC analysis, large database files are required but not included
in the package due to size constraints. The package includes automatic
download functionality.

### Automatic Download (Recommended)

``` bash
# Download SCENIC database files
./inst/scripts/download_scenic_files.sh

# Check if files are available
./inst/scripts/step4_scenic.sh --check-only
```

### Manual Download (Alternative)

If automatic download fails, manually download from Zenodo: 1.
`hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather`
(~1.2GB) 2. `motifs-v9-nr.hgnc-m0.001-o0.0.tbl` (~99MB)

Place downloaded files in `inst/extdata/` directory.

### Running Analysis with Auto-Download

``` bash
# Step 4 will automatically download dependencies
bash inst/scripts/step4_scenic.sh

# With custom options
bash inst/scripts/step4_scenic.sh --threads 20 --force
```

## Quick Start

``` r
library(HCCMultiOmics)

# Run prognostic model with ensemble ML
result <- run_prognostic_model(candidate_genes)

# Visualize results
plot_lasso_coef(result)
plot_rsf_vimp(result)
plot_xgb_imp(result)
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" alt="" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
