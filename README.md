<!-- README.md is generated from README.Rmd. Please edit that file -->

![](./man/figures/abstract_1.png)

# HCCMultiOmics

HCCMultiOmics provides an integrative analysis workflow for hepatocellular carcinoma (HCC) multi-omics data, spanning from bulk RNA-seq biomarker discovery to single-cell mechanism validation.

---

## Working Directory

All scripts must be executed from the `inst/scripts` directory:

```bash
cd /path/to/HCCMultiOmics/inst/scripts
```

Scripts rely on relative paths for input data and output files are saved to `hcc_output/`. Each step reads results from the previous step.

---

## Workflow Overview

The analysis workflow branches based on the cell type with highest target gene expression identified in Step 3:

```
Step 1: Bulk Analysis → Candidate Genes
         ↓
Step 2: Single-cell Visualization → Identify cell type with highest expression
         ↓
Step 3: Mechanism Analysis
         ↓
    ┌────┴────┐
    ↓         ↓
 Hepatocytes  Non-hepatocytes
 (Malignant)  (TME cells)
    ↓         ↓
Step 4-5    End
(SCENIC)    (Microenvironment
```

Workflow branches at Step 3:
- Hepatocytes (malignant cells): Continue to Step 4-5
- Non-hepatocytes (immune/stromal cells): Microenvironment analysis (workflow ends)

---

## Data Preparation

### Built-in Data

The following data are included in the package:
- TCGA expression matrix (`TCGA_Expr_Mat.rda`)
- TCGA clinical data (`TCGA_Clin_Data.rda`)
- TCGA LIHC differential genes (`TCGA_LIHC_DEGs.rda`)
- Built-in gene sets (Ferroptosis, Cuproptosis, etc.)

### Single-cell Data (Required for Step 2)

Option 1: Auto-download
```bash
Rscript step2_singlecell.R --gene YOUR_GENE --data auto
```

Option 2: Manual download from Zenodo: https://zenodo.org/records/18642056

Place `HCC_sc_data.rda` in the `inst/scripts/` directory.

### SCENIC Database Files (Required for Step 4)

Option 1: Auto-download
```bash
bash download_scenic_files.sh
```

Option 2: Manual download from Zenodo: https://zenodo.org/records/18901480

Place files in `inst/scripts/scenic_extdata/`:
- `hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather`
- `motifs-v9-nr.hgnc-m0.001-o0.0.tbl`
- `hs_hgnc_tfs.txt`

Python dependencies:
```bash
pip install pandas numpy pyscenic arboreto ctxcore
```

---

## Step-by-Step Guide

### Step 1: Bulk RNA-seq Analysis

**Purpose:** Build prognostic and diagnostic models to identify candidate genes associated with prognosis from the selected gene set.

**What it does:**
1. Intersects the selected gene set with TCGA LIHC differential genes
2. Trains three ensemble ML models (LASSO, RSF, XGBoost) for survival prediction
3. Identifies genes selected by all three models as final candidates
4. Builds a diagnostic model based on the candidate genes

**How to run:**

```bash
cd inst/scripts

# Built-in gene sets
Rscript step1_bulk_analysis.R --builtin 1

# Custom gene set
Rscript step1_bulk_analysis.R --custom genes.csv

# Set random seed
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

**Custom Gene Set Format (CSV):**

```csv
all_genes,drivers,suppressors,syml
GeneA,GeneB,GeneC,
GeneD,,GeneE,GeneF
```

All four columns required. Use empty values when no genes exist in a category.

**Output:** `hcc_output/step1/`

**Result Tables:**
- `candidate_genes.csv` - List of candidate genes
- `diag_result.rds` - Diagnostic model results

**Result Figures:**
- `venn_diagram.pdf` - Venn diagram of gene set intersection with TCGA DEGs
- `lasso_cv.pdf` - LASSO cross-validation plot
- `lasso_coef.pdf` - LASSO coefficient plot
- `rsf_process.pdf` - RSF survival curves
- `rsf_vimp.pdf` - RSF variable importance plot
- `xgb_importance.pdf` - XGBoost feature importance plot
- `model_venn.pdf` - Venn diagram of genes selected by 3 models
- `diag_heatmap.pdf` - Diagnostic model heatmap

**Next step:** Select one candidate gene for Step 2.

---

### Step 2: Single-cell Visualization

**Purpose:** Visualize target gene expression across cell types in single-cell data to determine the downstream analysis path.

**What it does:**
1. Generates diagnostic score bar plot from Step 1
2. Generates survival curve using TCGA data
3. Loads single-cell data and visualizes cell landscape (UMAP)
4. Visualizes target gene expression on UMAP
5. Calculates mean expression (Type 1) and positive cell ratio (Type 2) per cell type
6. Identifies cell type with highest expression

**How to run:**

```bash
cd inst/scripts

Rscript step2_singlecell.R --gene YOUR_GENE

# Auto-download single-cell data
Rscript step2_singlecell.R --gene YOUR_GENE --data auto
```

**Output:** `hcc_output/step2/`

**Result Tables:**
- `target_gene.rds` - Target gene name
- `top_cells.rds` - Cell types with highest expression

**Result Figures:**
- `sc_landscape.pdf` - Single-cell landscape
- `{GENE}_diag_score.pdf` - Diagnostic score by cell type
- `{GENE}_survival.pdf` - Survival curve
- `{GENE}_landscape.pdf` - Gene expression on UMAP
- `{GENE}_stat_type1.pdf` - Mean expression by cell type
- `{GENE}_stat_type2.pdf` - Positive cell ratio by cell type

**Next step:** Select Type 1 or Type 2 based on results. If highest expression in hepatocytes, proceed to Step 3.

---

### Step 3: Mechanism Analysis

**Purpose:** Analyze molecular mechanisms of the target gene. This step determines the downstream path based on cell type with highest gene expression.

**What it does:**
1. Reads target gene and top expressing cell type from Step 2
2. Branches based on cell type:

   **Hepatocytes path:**
   - Extracts hepatocytes and identifies malignant cells
   - Performs CytoTRACE stemness analysis
   - Performs trajectory/pseudotime analysis
   - Analyzes driver/suppressor/syml gene trends along pseudotime
   - Performs synthetic lethal gene screening
   - Generates SCENIC input for Step 4

   **Non-hepatocytes path:**
   - Performs ligand-receptor interaction analysis
   - Performs stromal cell function analysis
   - Performs immune microenvironment analysis

**How to run:**

```bash
cd inst/scripts

# Type 1: mean expression
Rscript step3_mechanism.R --type 1

# Type 2: positive cell ratio
Rscript step3_mechanism.R --type 2
```

**Output:** `hcc_output/step3/`

**Hepatocytes path - Result Figures:**
- `{GENE}_by_site.pdf` - Gene expression across sites
- `malignant_stats.pdf` - Malignant cell statistics
- `malignant_umap.pdf` - Malignant cell UMAP
- `malignant_clusters.pdf` - Malignant cell clustering
- `cyto_score.pdf` - CytoTRACE stemness score
- `cyto_violin.pdf` - Stemness score violin plot
- `trajectory.pdf` - Trajectory analysis
- `pseudotime.pdf` - Pseudotime distribution
- `driver_cor.pdf` - Driver gene correlations
- `driver_trend.pdf` - Driver gene trends
- `suppressor_cor.pdf` - Suppressor gene correlations
- `suppressor_trend.pdf` - Suppressor gene trends
- `syml_screening.pdf` - Synthetic lethal gene screening
- `{GENE}_{SYML}_density.pdf` - Density plot
- `tradeseq_heatmap.pdf` - TradeSeq heatmap
- `tradeseq_trends.pdf` - Gene expression trends
- `{GENE}_pseudotime_trend.pdf` - Target gene trend

**Result Tables:**
- `scenic_input.csv` - Expression matrix for SCENIC

**Next step:** Proceed to Step 4.

---

**Non-hepatocytes path - Result Figures:**
- `ligand_bubble.pdf` - Ligand-receptor interaction
- `stromal_function.pdf` - Stromal cell function
- `immune_landscape.pdf` - Immune cell landscape
- `immune_scatter_grid.pdf` - Immune cell scatter plots

Workflow ends here.

---

### Step 4: SCENIC Analysis

**Purpose:** Build gene regulatory network using SCENIC to identify transcription factors regulating the target gene.

**What it does:**
1. Constructs co-expression network using GRNBoost2
2. Prunes network with motif enrichment (CTX)
3. Identifies TF regulons
4. Calculates AUC scores for each regulon

**How to run:**

```bash
cd inst/scripts

bash download_scenic_files.sh
bash step4_scenic.sh 10
```

Python environment with pyscenic is required.

**Output:** `hcc_output/step4/`

**Result Tables:**
- `regulons.csv` - TF regulons
- `AUC_matrix.csv` - AUC scores

**Next step:** Proceed to Step 5.

---

### Step 5: SCENIC Downstream Analysis

**Purpose:** Identify key transcription factors and analyze regulatory relationships between TFs, target gene, and downstream gene.

**What it does:**
1. Loads SCENIC AUC matrix and regulons from Step 4
2. Performs TF enrichment analysis
3. Screens for mediator TFs
4. Visualizes TF-target relationships

**How to run:**

```bash
cd inst/scripts

Rscript step5_scenic_analysis.R --target YOUR_GENE --downstream DOWNSTREAM_GENE
```

**Output:** `hcc_output/step5/`

**Result Figures:**
- TF enrichment heatmaps
- Mediator TF analysis plots
- Network visualizations

Workflow complete.

---

## Quick Reference

| Step | Command | Input | Output |
|------|---------|-------|--------|
| 1 | `Rscript step1_bulk_analysis.R --builtin 1` | Gene set | candidate_genes.csv |
| 2 | `Rscript step2_singlecell.R --gene GENE` | Single-cell data | top_cells.rds |
| 3 | `Rscript step3_mechanism.R --type 1/2` | target_gene.rds | scenic_input.csv |
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
