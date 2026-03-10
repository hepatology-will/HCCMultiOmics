#!/bin/bash
# Step 4: SCENIC Analysis
# Runs SCENIC analysis using local Python script

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# Defaults
THREADS=10
INPUT="inst/scripts/hcc_output/step3/scenic_input.csv"
OUTPUT="inst/scripts/hcc_output/step4"
DATA="scenic_extdata"

# Help
if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    echo "Usage: $0 [threads]"
    echo "  threads: Number of CPU threads (default: 10)"
    echo ""
    echo "Example: $0 20"
    exit 0
fi

if [ -n "$1" ]; then
    THREADS="$1"
fi

echo "=== Step 4: SCENIC Analysis ==="
echo "Threads: $THREADS"
echo ""

# Check input file
echo "1. Checking input file..."
if [ ! -f "$INPUT" ]; then
    echo "✗ Input file not found: $INPUT"
    echo "  Run step3 first to generate scenic_input.csv"
    exit 1
fi
INPUT_SIZE=$(du -h "$INPUT" | cut -f1)
echo "✓ Input file: $INPUT"
echo "  Size: $INPUT_SIZE"
echo ""

# Check SCENIC database files
echo "2. Checking SCENIC database files..."
FILES=(
    "$DATA/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
    "$DATA/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
    "$DATA/hs_hgnc_tfs.txt"
)

for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✓ $(basename "$file")"
    else
        echo "  ✗ $(basename "$file")"
        echo "  Run ./download_scenic_files.sh first"
        exit 1
    fi
done

echo ""
echo "✓ All SCENIC database files found"
echo ""

# Create output directory
mkdir -p "$OUTPUT"

# Prepare Python script with correct paths
echo "3. Preparing SCENIC analysis..."
PYTHON_SCRIPT="$OUTPUT/run_scenic_fixed.py"

# Get absolute paths for Python script
ABS_INPUT="$(readlink -f "$INPUT" || echo "$INPUT")"
ABS_DATA="$(readlink -f "$DATA" || echo "$DATA")"
ABS_OUTPUT="$(readlink -f "$OUTPUT" || echo "$OUTPUT")"

# Create fixed Python script
cat > "$PYTHON_SCRIPT" << EOF
import os
import pandas as pd
import numpy as np
from arboreto.algo import grnboost2

# --- 修正点 1: 新版 API 路径 ---
try:
    from ctxcore.rnkdb import FeatherRankingDatabase 
    from pyscenic.utils import modules_from_adjacencies, load_motifs
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.aucell import aucell
except ImportError as e:
    print(f"Error importing pyscenic modules: {e}")
    print("Make sure pyscenic is installed: pip install pyscenic")
    exit(1)

# ================= 配置区域 =================
F_EXPRESSION = "$ABS_INPUT"
F_TFS = "$ABS_DATA/hs_hgnc_tfs.txt"
F_MOTIFS = "$ABS_DATA/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
F_DB = "$ABS_DATA/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"

F_OUTPUT_CSV = "$ABS_OUTPUT/AUC_matrix.csv"
F_OUTPUT_REGULONS = "$ABS_OUTPUT/regulons.csv"

N_CPU = $THREADS
# ===========================================

def load_tf_names(tf_file):
    with open(tf_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def run():
    print("Step 1: Loading Expression Matrix...")
    ex_matrix = pd.read_csv(F_EXPRESSION, index_col=0)
    print(f"Matrix shape: {ex_matrix.shape} (Cells x Genes)")
    
    tf_names = load_tf_names(F_TFS)
    print(f"Loaded {len(tf_names)} TFs.")

    print("Step 2: Inferring Co-expression Network (GRNBoost2)...")
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True, num_workers=N_CPU)

    print("Step 3: Pruning Network with Motifs (CTX)...")
    modules = modules_from_adjacencies(adjacencies, ex_matrix)
    
    dbs = [FeatherRankingDatabase(fname=F_DB, name="hg38_10kb")]
    
    df = prune2df(dbs, modules, F_MOTIFS)
    regulons = df2regulons(df)
    print(f"Found {len(regulons)} regulons.")
    
    # Save regulons
    df.to_csv(F_OUTPUT_REGULONS, index=False)
    print(f"Saved regulons to: {F_OUTPUT_REGULONS}")

    print("Step 4: Calculating AUC Scores...")
    auc_mtx = aucell(ex_matrix, regulons, num_workers=N_CPU)
    
    # Save AUC matrix
    auc_mtx.to_csv(F_OUTPUT_CSV)
    print(f"Saved AUC matrix to: {F_OUTPUT_CSV}")
    
    print("\\n✓ SCENIC analysis completed successfully!")
    print(f"  - Regulons: {F_OUTPUT_REGULONS}")
    print(f"  - AUC matrix: {F_OUTPUT_CSV}")

if __name__ == "__main__":
    run()
EOF

echo "4. Checking Python environment..."
if ! python3 -c "import pandas, numpy, pyscenic" 2>/dev/null; then
    echo "✗ Python environment not ready for SCENIC"
    echo ""
    echo "Required packages:"
    echo "  pandas, numpy, pyscenic, arboreto, ctxcore"
    echo ""
    echo "Install with:"
    echo "  pip install pandas numpy pyscenic arboreto ctxcore"
    echo "  or"
    echo "  conda install -c conda-forge -c bioconda pyscenic"
    exit 1
fi

echo "✓ Python environment ready"
echo ""
echo "5. Running SCENIC analysis..."
echo "   This may take several hours..."
echo ""

# Run Python script
if python3 "$PYTHON_SCRIPT"; then
    echo ""
    echo "✓ SCENIC analysis completed!"
    echo "  Results saved to: $OUTPUT"
    echo ""
    echo "Output files:"
    ls -lh "$OUTPUT"/*.csv 2>/dev/null || echo "  (No CSV files found)"
else
    echo ""
    echo "✗ SCENIC analysis failed"
    echo "  Check Python environment and dependencies"
    exit 1
fi

echo ""
echo "========== Step 4 Complete! =========="
echo "Next step: Run step5_scenic_analysis.R"
echo "  Rscript step5_scenic_analysis.R --target EZH2 --downstream SLC7A11"