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
1. Intersects your selected gene set with TCGA LIHC differential genes to find candidate genes
2. Trains three ensemble ML models (LASSO, RSF, XGBoost) for survival prediction
3. Identifies genes selected by all three models as final candidate genes
4. Builds a diagnostic model based on the candidate genes

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

**Custom Gene Set Format (CSV only):**

```csv
all_genes,drivers,suppressors,syml
GeneA,GeneB,GeneC,
GeneD,,GeneE,GeneF
```

All four columns required. Use empty values if no genes in that category.

**Output:** `hcc_output/step1/`

**Result Tables:**
- `candidate_genes.csv` - List of candidate genes for next step
- `diag_result.rds` - Diagnostic model results (RDS format)

**Result Figures:**
- `venn_diagram.pdf` - Venn diagram showing gene set intersection with TCGA DEGs
- `lasso_cv.pdf` - LASSO cross-validation plot showing optimal lambda
- `lasso_coef.pdf` - LASSO coefficient plot showing selected genes
- `rsf_process.pdf` - Random Survival Forest survival curves
- `rsf_vimp.pdf` - RSF variable importance plot
- `xgb_importance.pdf` - XGBoost feature importance plot
- `model_venn.pdf` - Venn diagram of genes selected by all 3 models
- `diag_heatmap.pdf` - Diagnostic model heatmap

**How to interpret:**
1. Check `venn_diagram.pdf` to see overlap between your gene set and TCGA DEGs
2. Review model performance plots (`lasso_cv.pdf`, `rsf_process.pdf`) to assess model quality
3. Select candidate genes from `candidate_genes.csv` based on model selection
4. Choose one target gene for Step 2

**Next step:** Select one candidate gene from the output and proceed to Step 2.

---

### Step 2: Single-cell Visualization

**Purpose:** Visualize target gene expression across different cell types in single-cell data to determine the downstream analysis path.

**What it does:**
1. Generates diagnostic score bar plot from Step 1 diagnostic model results
2. Generates survival curve using TCGA data to show prognostic value
3. Loads single-cell data and visualizes overall cell landscape (UMAP)
4. Visualizes target gene expression on UMAP (FeaturePlot + Violin)
5. Calculates mean expression (Type 1) and positive cell ratio (Type 2) per cell type
6. Identifies which cell type has the highest expression

**How to run:**

```bash
cd inst/scripts

# Replace YOUR_GENE with your selected gene name from Step 1
Rscript step2_singlecell.R --gene YOUR_GENE

# Auto-download single-cell data from Zenodo if not available locally
Rscript step2_singlecell.R --gene YOUR_GENE --data auto
```

**Output:** `hcc_output/step2/`

**Result Tables:**
- `target_gene.rds` - Target gene name (for Step 3)
- `top_cells.rds` - Cell types with highest expression (for Step 3)

**Result Figures:**
- `sc_landscape.pdf` - Overall single-cell landscape showing all cell types
- `{GENE}_diag_score.pdf` - Diagnostic score comparison across cell types
- `{GENE}_survival.pdf` - Survival curve stratified by gene expression
- `{GENE}_landscape.pdf` - Gene expression on UMAP (FeaturePlot + Violin)
- `{GENE}_stat_type1.pdf` - Mean expression by cell type (Type 1)
- `{GENE}_stat_type2.pdf` - Positive cell ratio by cell type (Type 2)

**How to interpret:**
1. Review `sc_landscape.pdf` to understand cell type composition
2. Check `{GENE}_stat_type1.pdf` and `{GENE}_stat_type2.pdf` to identify which cell type has highest expression
3. Look at `{GENE}_survival.pdf` for prognostic value
4. Choose the type (1 or 2) that shows better separation between cell types
5. If highest expression is in **hepatocytes** (malignant cells) → Proceed to Step 3
6. If highest expression is in **non-hepatocytes** (immune/stromal cells) → Microenvironment analysis ends here

**Key Decision:** Based on which cell type shows highest expression:
- **Hepatocytes** → Step 3 (Mechanism + SCENIC workflow)
- **Non-hepatocytes** → Workflow ends (microenvironment analysis)

---

### Step 3: Mechanism Analysis

**Purpose:** Analyze the molecular mechanism of your target gene. This step automatically determines the downstream path based on which cell type has the highest gene expression (from Step 2).

**What it does:**
1. Automatically reads target gene and top expressing cell type from Step 2 results
2. Determines workflow based on cell type:
   - **If Hepatocytes (malignant cells):** Workflow 1
     - Extracts hepatocytes and identifies malignant cells
     - Performs CytoTRACE stemness analysis
     - Performs trajectory/pseudotime analysis
     - Analyzes driver/suppressor/syml gene trends along pseudotime
     - Performs synthetic lethal gene screening
     - Generates SCENIC input for Step 4
   - **If Non-hepatocytes (TME cells):** Workflow 2
     - Performs ligand-receptor interaction analysis
     - Performs stromal cell function analysis
     - Performs immune microenvironment analysis

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

This step has two possible workflows based on cell type:

#### Workflow 1: Hepatocytes (Malignant Cells)

**Result Figures:**
- `{GENE}_by_site.pdf` - Gene expression across different sites
- `malignant_stats.pdf` - Malignant cell statistics
- `malignant_umap.pdf` - Malignant cell UMAP visualization
- `malignant_clusters.pdf` - Malignant cell clustering
- `cyto_score.pdf` - CytoTRACE stemness score
- `cyto_violin.pdf` - Stemness score violin plot
- `trajectory.pdf` - Trajectory analysis (pseudotime)
- `pseudotime.pdf` - Pseudotime distribution
- `driver_cor.pdf` - Driver gene correlations
- `driver_trend.pdf` - Driver gene trends along pseudotime
- `suppressor_cor.pdf` - Suppressor gene correlations
- `suppressor_trend.pdf` - Suppressor gene trends along pseudotime
- `syml_screening.pdf` - Synthetic lethal gene screening
- `{GENE}_{SYML}_density.pdf` - Density plot for top synthetic lethal gene
- `tradeseq_heatmap.pdf` - TradeSeq heatmap
- `tradeseq_trends.pdf` - Gene expression trends
- `{GENE}_pseudotime_trend.pdf` - Target gene trend along pseudotime

**Result Tables:**
- `scenic_input.csv` - Expression matrix for SCENIC analysis

**How to interpret:**
1. Review `malignant_umap.pdf` and `malignant_clusters.pdf` to understand malignant cell populations
2. Check `cyto_score.pdf` and `cyto_violin.pdf` for stemness analysis
3. Review `trajectory.pdf` and `pseudotime.pdf` for developmental trajectory
4. Analyze driver/suppressor gene trends along pseudotime
5. Use `scenic_input.csv` for Step 4 SCENIC analysis

**Next step:** Proceed to Step 4 (SCENIC Analysis)

---

#### Workflow 2: Non-hepatocytes (TME Cells)

**Result Figures:**
- `ligand_bubble.pdf` - Ligand-receptor interaction bubble plot
- `stromal_function.pdf` - Stromal cell function analysis
- `immune_landscape.pdf` - Immune cell landscape visualization
- `immune_scatter_grid.pdf` - Immune cell scatter plots

**How to interpret:**
1. Review `ligand_bubble.pdf` for ligand-receptor interactions
2. Check `stromal_function.pdf` for stromal cell functional analysis
3. Review `immune_landscape.pdf` and `immune_scatter_grid.pdf` for immune microenvironment

**Workflow ends here** - This path analyzes tumor microenvironment interactions.

---

### Step 4: SCENIC Analysis

**Purpose:** Build gene regulatory network using SCENIC to identify transcription factors (TFs) that may regulate your target gene.

**What it does:**
1. Constructs co-expression network using GRNBoost2 algorithm
2. Prunes network with motif enrichment using CTX (ctxcore)
3. Identifies TF regulons (TF-target gene pairs)
4. Calculates AUC (Area Under the Curve) scores for each regulon across single cells

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

**Result Tables:**
- `regulons.csv` - Transcription factor regulons (TF-target gene relationships)
- `AUC_matrix.csv` - AUC scores for each regulon in each cell

**How to interpret:**
1. `regulons.csv` contains TF-target gene pairs identified by SCENIC
2. `AUC_matrix.csv` shows regulon activity scores across single cells
3. These files are used for downstream TF enrichment analysis in Step 5

**Next step:** Proceed to Step 5 for downstream analysis.

---

### Step 5: SCENIC Downstream Analysis

**Purpose:** Identify key transcription factors and analyze regulatory relationships between TFs, your target gene, and a downstream gene of interest.

**What it does:**
1. Loads SCENIC AUC matrix and regulons from Step 4
2. Performs TF enrichment analysis to find TFs with high activity
3. Screens for mediator TFs that connect your target gene to the downstream gene
4. Visualizes TF-target gene regulatory relationships as heatmaps and network plots

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

**Result Figures:**
- TF enrichment heatmaps showing TF activity
- Mediator TF analysis plots
- Network visualizations showing TF-target gene regulatory relationships

**How to interpret:**
1. Review TF enrichment heatmaps to identify TFs with high activity
2. Check mediator TF analysis for intermediate TFs between your target and downstream gene
3. Network plots show regulatory relationships between TFs and target genes

**Workflow complete!**

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
