#!/usr/bin/env bash
# Executes prodigal for all folder containing genomes grouped by clade

set -euo pipefail

for dir in group_genomes_*; do

    [[ -d "$dir" ]] || continue

    # Extract suffix
    suffix="${dir#group_genomes_}"

    # Directory output name
    outdir="group_prodigal_genomes_${suffix}"

    echo "Running prodigal on $dir → $outdir"
    bash /home/estudiante2/JMF/final_analysis/prodigal/prodigal_analysis.sh "$dir" "$outdir"
done
