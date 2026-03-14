<!-- README.md is generated from README.Rmd. Please edit that file -->

![](./man/figures/abstract_1.png)

# HCCMultiOmics

[![R-CMD-check](https://github.com/hepatology-will/HCCMultiOmics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hepatology-will/HCCMultiOmics/actions/workflows/R-CMD-check.yaml)

HCCMultiOmics is an integrative analysis framework for multiscale biological discovery in Hepatocellular Carcinoma (HCC). It provides a complete workflow from bulk RNA-seq biomarker discovery to single-cell mechanism validation.

## Installation

```r
# Install devtools if not already installed
install.packages("devtools")

# Install HCCMultiOmics
devtools::install_github("hepatology-will/HCCMultiOmics")
```

## Workflow Overview

The package provides a **5-step workflow**:

```
Step 1: Bulk Analysis (Machine Learning)
    ↓
Step 2: Single-cell Gene Expression Visualization  
    ↓
Step 3: Mechanism Analysis
    ↓
Step 4: SCENIC Analysis (Python)
    ↓
Step 5: SCENIC Downstream Analysis
```

## Step-by-Step Guide

### Step 1: Bulk RNA-seq Analysis

Select candidate genes using ensemble machine learning (LASSO, RSF, XGBoost).

```bash
# Using built-in gene sets
Rscript inst/scripts/step1_bulk_analysis.R --builtin 1

# Using custom gene set
Rscript inst/scripts/step1_bulk_analysis.R --custom genes.csv
```

**Built-in Gene Sets:**
- 1. Ferroptosis_FerrDb
- 2. Cuproptosis_FerrDb
- 3. Disulfidptosis
- 4. Autosis
- 5. immunogonic_cell_death
- 6. mitotic_death
- 7. parthanatos

**Output:** Candidate gene list in `hcc_output/step1/`

---

### Step 2: Single-cell Visualization

Visualize target gene expression in single-cell data.

```bash
# Using local single-cell data
Rscript inst/scripts/step2_singlecell.R --gene YOUR_GENE

# Auto-download single-cell data from Zenodo
Rscript inst/scripts/step2_singlecell.R --gene YOUR_GENE --data auto
```

**Output:** Expression plots in `hcc_output/step2/`

You will get two analysis types:
- **Type 1**: Mean expression
- **Type 2**: Positive cell ratio

Choose the type with better separation for downstream analysis.

---

### Step 3: Mechanism Analysis

Analyze gene mechanisms in malignant hepatocytes.

```bash
# Using mean expression (type 1)
Rscript inst/scripts/step3_mechanism.R --type 1

# Using positive cell ratio (type 2)
Rscript inst/scripts/step3_mechanism.R --type 2
```

**Output:** Mechanism analysis results in `hcc_output/step3/`

---

### Step 4: SCENIC Analysis

Run SCENIC for transcription factor regulatory network analysis.

```bash
# First, download SCENIC database files (required)
./inst/scripts/download_scenic_files.sh

# Run SCENIC (requires Python environment with pyscenic)
bash inst/scripts/step4_scenic.sh 10
```

**Python Requirements:**
```bash
pip install pandas numpy pyscenic arboreto ctxcore
```

**Output:** AUC matrix and regulons in `hcc_output/step4/`

---

### Step 5: SCENIC Downstream Analysis

Analyze TF enrichment and visualize results.

```bash
Rscript inst/scripts/step5_scenic_analysis.R --target EZH2 --downstream SLC7A11
```

**Output:** TF enrichment plots in `hcc_output/step5/`

---

## Data Requirements

### Required Data Files

Place the following files in the working directory:

| File | Description |
|------|-------------|
| `TCGA_Expr_Mat.rda` | TCGA expression matrix |
| `TCGA_Clin_Data.rda` | TCGA clinical data |
| `TCGA_LIHC_DEGs.rda` | TCGA LIHC differential genes |
| `HCC_sc_data.rda` | Single-cell data (or auto-download) |

### SCENIC Database Files

Download via:
```bash
./inst/scripts/download_scenic_files.sh
```

Or manually from [Zenodo](https://zenodo.org/records/18642056), place in `inst/scripts/scenic_extdata/`.

## Package Functions

Key functions for programmatic use:

```r
# Run prognostic model
result <- run_prognostic_model(candidate_genes, seed = 123)

# Run diagnostic model
diag_result <- run_diagnostic_model(candidate_genes)

# Visualization
plot_lasso_cv(result)
plot_lasso_coef(result)
plot_rsf_process(result)
plot_rsf_vimp(result)
plot_xgb_imp(result)
plot_gene_survival(gene = "EZH2")
plot_sc_gene(seurat_object, gene = "EZH2")
```

## Citation

```bibtex
@software{HCCMultiOmics,
  author = {Zhuo Chen},
  title = {HCCMultiOmics: An Integrative Analysis Framework for HCC Multi-omics},
  year = {2025},
  url = {https://github.com/hepatology-will/HCCMultiOmics}
}
```

## License

MIT License
