#!/bin/bash
set -euo pipefail

BAM_DIR="bam_subset"
OUTPUT_TSV="downsample_input.tsv"
COVERAGES="[1,2,5,10,15]"

echo -e "pair_id\tpath\tdownsample_cov" > "$OUTPUT_TSV"

# Loop over .bam files only (exclude .bai)
for bamfile in "$BAM_DIR"/*.bam; do
    if [[ "$bamfile" == *.bam.bai ]]; then
        continue
    fi
    filename=$(basename "$bamfile")
    pair_id="${filename%.bam}"

    echo -e "${pair_id}\t${bamfile}\t${COVERAGES}" >> "$OUTPUT_TSV"
done

echo "✅ Generated $OUTPUT_TSV"