<!-- README.md is generated from README.Rmd. Please edit that file -->



# HCCMultiOmics

<!-- badges: start -->
[![R-CMD-check](https://github.com/hepatology-will/HCCMultiOmics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hepatology-will/HCCMultiOmics/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

HCCMultiOmics is an integrative AutoML framework for multiscale biological discovery in Hepatocellular Carcinoma (HCC). The package provides a comprehensive and automated workflow designed to translate multi-omics data into actionable biological mechanisms.

## Features

- **Module A: Bulk Analysis**: Ensemble machine learning engine (LASSO, RSF, XGBoost) for robust prognostic biomarker discovery
- **Module B: Single-cell Analysis**: Adaptive single-cell mechanistic decoding module to profile intracellular trajectories, stemness, and intercellular communication
- **Integrated Workflow**: Seamless transition from bulk to single-cell analysis
- **Automated Visualization**: Comprehensive plotting functions for results interpretation

## Installation

You can install the development version of HCCMultiOmics from [GitHub](https://github.com/hepatology-will/HCCMultiOmics) with:

``` r
# install.packages("devtools")
devtools::install_github("hepatology-will/HCCMultiOmics")
```

## Data Requirements

For SCENIC analysis, large database files are required but not included in the package due to size constraints. The package includes automatic download functionality.

### Automatic Download (Recommended)

```bash
# Download SCENIC database files
./inst/scripts/download_scenic_files.sh

# Check if files are available
./inst/scripts/step4_scenic.sh --check-only
```

### Manual Download (Alternative)

If automatic download fails, manually download from Zenodo:
1. `hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather` (~1.2GB)
2. `motifs-v9-nr.hgnc-m0.001-o0.0.tbl` (~99MB)

Place downloaded files in `inst/extdata/` directory.

### Running Analysis with Auto-Download

```bash
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

## Documentation

For detailed documentation and examples, please visit the [package website](https://github.com/hepatology-will/HCCMultiOmics).

## Citation

If you use HCCMultiOmics in your research, please cite:

``` bibtex
@software{HCCMultiOmics,
  author = {Zhuo Chen},
  title = {HCCMultiOmics: An Integrative AutoML Framework for Multiscale Biological Discovery},
  year = {2025},
  url = {https://github.com/hepatology-will/HCCMultiOmics}
}
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This package is licensed under the MIT License. See the LICENSE file for details.
