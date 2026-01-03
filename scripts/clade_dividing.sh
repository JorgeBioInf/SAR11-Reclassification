#!/usr/bin/env bash
# Divides test_genomes/ directory into subdirectories with the genomes from each clade

set -euo pipefail

# Check argument
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <input_directory>"
    exit 1
fi

INPUT_DIR="$1"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: '$INPUT_DIR' is not a directory"
    exit 1
fi

for f in "$INPUT_DIR"/*.fa; do

    [[ -e "$f" ]] || continue

    filename=$(basename "$f")

    # Extract clade
    clade="${filename##*_}"
    clade="${clade%.fa}"
    out_dir="clade_$clade"

    mkdir -p "$INPUT_DIR/$out_dir"
    mv "$f" "$INPUT_DIR/$out_dir/"
done