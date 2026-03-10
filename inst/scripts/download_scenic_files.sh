#!/bin/bash
# Simple script to download SCENIC database files to scenic_extdata/ in working directory

set -e

DEST_DIR="scenic_extdata"
FILES=(
    "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
    "motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
    "hs_hgnc_tfs.txt"
)

URLS=(
    "https://zenodo.org/records/18901480/files/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
    "https://zenodo.org/records/18901480/files/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
    "https://zenodo.org/records/18901480/files/hs_hgnc_tfs.txt"
)

echo "Setting up SCENIC database files in $DEST_DIR/"

# Create directory if it doesn't exist
if [ ! -d "$DEST_DIR" ]; then
    mkdir -p "$DEST_DIR"
    echo "Created directory: $DEST_DIR"
fi

# Download files
for i in "${!FILES[@]}"; do
    FILE="${FILES[$i]}"
    URL="${URLS[$i]}"
    DEST="$DEST_DIR/$FILE"
    
    if [ -f "$DEST" ]; then
        echo "File already exists: $FILE"
    else
        echo "Downloading: $FILE"
        wget -q --show-progress -O "$DEST" "$URL"
        
        if [ $? -eq 0 ]; then
            echo "  ✓ Downloaded: $FILE"
        else
            echo "  ✗ Failed to download: $FILE"
            exit 1
        fi
    fi
done

echo ""
echo "All SCENIC database files are ready in $DEST_DIR/"
echo "Files:"
ls -lh "$DEST_DIR/"