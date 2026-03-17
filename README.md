<!-- README.md is generated from README.Rmd. Please edit that file -->

![](./man/figures/abstract_1.png)

# HCCMultiOmics

HCCMultiOmics provides an integrative analysis workflow for hepatocellular carcinoma (HCC) multi-omics data, spanning from bulk RNA-seq biomarker discovery to single-cell mechanism validation.

---

## Important: Working Directory

**All scripts must be run from the `inst/scripts` directory:**

```bash
cd /path/to/HCCMultiOmics/inst/scripts
```

This is critical because:
- Scripts read input data from relative paths
- Output files are saved to `hcc_output/`
- Step N reads results from Step N-1

---

## Workflow Overview

The analysis workflow is **decision-tree based** - the downstream path is determined by Step 3 based on the cell type with highest target gene expression:

```
Step 1: Bulk Analysis → Candidate Genes
         ↓
Step 2: Single-cell Visualization → Identify cell type with highest expression
         ↓
Step 3: Mechanism Analysis (determined by cell type)
         ↓
    ┌────┴────┐
    ↓         ↓
 Hepatocytes  Non-hepatocytes
 (Malignant)  (TME cells)
    ↓         ↓
Step 4-5    End
(SCENIC)    (Microenvironment
 +Downstream  analysis)
```

**Key Decision (Step 3):**
- If target gene is highly expressed in **malignant hepatocytes** → Continue to Step 4-5 (SCENIC analysis)
- If target gene is highly expressed in **non-hepatocytes** (immune/stromal cells) → Workflow ends (microenvironment analysis)

---

## Data Preparation

### Built-in Data (No Preparation Needed)

The following data are **included in the package**:
- TCGA expression matrix (`TCGA_Expr_Mat.rda`)
- TCGA clinical data (`TCGA_Clin_Data.rda`)
- TCGA LIHC differential genes (`TCGA_LIHC_DEGs.rda`)
- Built-in gene sets (Ferroptosis, Cuproptosis, etc.)

### Single-cell Data (Required for Step 2)

**Option 1: Auto-download (Recommended)**
```bash
Rscript step2_singlecell.R --gene YOUR_GENE --data auto
```

**Option 2: Manual Download**
Download from Zenodo: https://zenodo.org/records/18642056

Place `HCC_sc_data.rda` in the `inst/scripts/` directory.

### SCENIC Database Files (Required for Step 4)

**Option 1: Auto-download (Recommended)**
```bash
bash download_scenic_files.sh
```

**Option 2: Manual Download**
Download from [Zenodo](https://zenodo.org/records/18901480) and place in `inst/scripts/scenic_extdata/`:
- `hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather`
- `motifs-v9-nr.hgnc-m0.001-o0.0.tbl`
- `hs_hgnc_tfs.txt`

Requires Python environment:
```bash
pip install pandas numpy pyscenic arboreto ctxcore
```

---

## Step-by-Step Guide

### Step 1: Bulk RNA-seq Analysis

**Purpose:** Identify candidate prognostic genes using machine learning on TCGA bulk RNA-seq data.

**What it does:**
1. Intersects your selected gene set with TCGA differential genes
2. Trains three ensemble ML models (LASSO, RSF, XGBoost)
3. Identifies candidate genes that are prognostic
4. Builds diagnostic model for the candidate genes

**How to run:**

```bash
cd inst/scripts

# Using built-in gene sets
Rscript step1_bulk_analysis.R --builtin 1

# Using custom gene set (see format below)
Rscript step1_bulk_analysis.R --custom genes.csv

# Set random seed for reproducibility
Rscript step1_bulk_analysis.R --builtin 1 --seed 123
```

**Built-in Gene Sets:**
| Index | Gene Set |
|-------|----------|
| 1 | Ferroptosis_FerrDb |
| 2 | Cuproptosis_FerrDb |
| 3 | Disulfidptosis |
| 4 | Autosis |
| 5 | immunogonic_cell_death |
| 6 | mitotic_death |
| 7 | parthanatos |

**Custom Gene Set Format:**

1. CSV/TXT file with columns (all required):

```csv
all_genes,drivers,suppressors,syml
GeneA,GeneB,GeneC,
GeneD,,GeneE,GeneF
```

All columns required. Use empty values if no genes in that category.

2. RDS/RData file: Save GeneSet object from the package

3. Plain text (one gene per line):
```
GeneA
GeneB
GeneC
```

**Output:** `hcc_output/step1/`
- `candidate_genes.csv` - List of candidate genes for next step
- `diag_result.rds` - Diagnostic model results
- Various PDF plots (model performance, Venn diagrams)

**Next step:** Select one candidate gene from the output and proceed to Step 2.

---

### Step 2: Single-cell Visualization

**Purpose:** Visualize target gene expression across different cell types in single-cell data to determine the downstream analysis path.

**What it does:**
1. Shows target gene expression on UMAP (FeaturePlot + Violin)
2. Calculates mean expression (Type 1) and positive cell ratio (Type 2) per cell type
3. Identifies which cell type has the highest expression

**How to run:**

```bash
cd inst/scripts

# Replace YOUR_GENE with your selected gene name from Step 1
Rscript step2_singlecell.R --gene YOUR_GENE

# Auto-download single-cell data from Zenodo if not available locally
Rscript step2_singlecell.R --gene YOUR_GENE --data auto
```

**Output:** `hcc_output/step2/`
- `sc_landscape.pdf` - Overall single-cell landscape
- `{GENE}_landscape.pdf` - Gene expression on UMAP
- `{GENE}_stat_type1.pdf` - Mean expression by cell type
- `{GENE}_stat_type2.pdf` - Positive cell ratio by cell type
- `top_cells.rds` - Cell types with highest expression (for Step 3)
- `target_gene.rds` - Target gene name (for Step 3)

**Key Decision:** Look at the output plots to decide which type (1 or 2) shows better separation, then proceed to Step 3.

---

### Step 3: Mechanism Analysis

**Purpose:** Analyze the molecular mechanism of your target gene. This step determines the downstream path based on cell type analysis.

**What it does:**
1. Analyzes trajectory, stemness, and intercellular communication
2. Determines if target gene is in hepatocytes or non-hepatocytes
3. If in hepatocytes: generates input for SCENIC analysis (Step 4)
4. If in non-hepatocytes: outputs microenvironment analysis results (workflow ends)

**How to run:**

```bash
cd inst/scripts

# Use Type 1 if mean expression showed better separation
Rscript step3_mechanism.R --type 1

# Use Type 2 if positive cell ratio showed better separation
Rscript step3_mechanism.R --type 2
```

**Important:** This step automatically reads:
- Target gene from `hcc_output/step2/target_gene.rds`
- Top expressing cell types from `hcc_output/step2/top_cells.rds`

**Output:** `hcc_output/step3/`
- Trajectory analysis plots
- Stemness analysis results
- Ligand-receptor communication analysis
- `scenic_input.csv` - Expression matrix for SCENIC (only if in hepatocytes path)

**Next step:**
- If target gene in hepatocytes → Proceed to Step 4
- If target gene in non-hepatocytes → Workflow ends here

---

### Step 4: SCENIC Analysis

**Purpose:** Build gene regulatory network using SCENIC to identify transcription factors (TFs) regulating your target gene.

**What it does:**
1. Constructs co-expression network (GRNBoost2)
2. Prunes network with motif enrichment (CTX)
3. Calculates AUC scores for regulons

**How to run:**

```bash
cd inst/scripts

# First ensure SCENIC database files are downloaded
bash download_scenic_files.sh

# Run SCENIC (10 threads, adjust as needed)
bash step4_scenic.sh 10
```

**Important:** Requires Python environment with pyscenic installed.

**Output:** `hcc_output/step4/`
- `regulons.csv` - TF regulons
- `AUC_matrix.csv` - AUC scores for each cell

**Next step:** Proceed to Step 5 for downstream analysis.

---

### Step 5: SCENIC Downstream Analysis

**Purpose:** Identify key transcription factors and analyze regulatory relationships.

**What it does:**
1. Performs TF enrichment analysis
2. Screens for mediator TFs
3. Visualizes TF-target gene relationships
4. Generates heatmaps and network plots

**How to run:**

```bash
cd inst/scripts

# --target: Your target gene from Step 2
# --downstream: Downstream gene of interest
Rscript step5_scenic_analysis.R --target YOUR_GENE --downstream DOWNSTREAM_GENE
```

**Example:**
```bash
Rscript step5_scenic_analysis.R --target EZH2 --downstream SLC7A11
```

**Output:** `hcc_output/step5/`
- TF enrichment heatmaps
- Mediator TF analysis plots
- Network visualizations

---

## Quick Reference

| Step | Command | Input | Output |
|------|---------|-------|--------|
| 1 | `Rscript step1_bulk_analysis.R --builtin 1` | Gene set (built-in or custom) | candidate_genes.csv |
| 2 | `Rscript step2_singlecell.R --gene GENE` | Single-cell data | top_cells.rds |
| 3 | `Rscript step3_mechanism.R --type 1/2` | target_gene.rds | scenic_input.csv (if hepatocytes) |
| 4 | `bash step4_scenic.sh 10` | scenic_input.csv | regulons.csv |
| 5 | `Rscript step5_scenic_analysis.R --target G --downstream G2` | AUC_matrix.csv | TF analysis |

---

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
