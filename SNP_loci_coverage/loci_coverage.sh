#!/bin/bash
set -euo pipefail

# Input files and directories
LOCI_FILE="loci.tsv"
BAM_DIR="downsampled_bams/"
OUTPUT="locus_depths.tsv"

# Loci to analyze
declare -A LOCI=(
    ["chr1:11796321"]=1
    ["chr19:44908684"]=1
    ["chr19:44908822"]=1
)

# Initialize output
echo -e "sample\tchrom\tpos\tdepth" > "$OUTPUT"

# First, verify we can find chr19 files
echo "=== VERIFYING FILES ==="
ls -la "$BAM_DIR"/*chr19*.bam || echo "No chr19 files found with basic pattern"

# Alternative method to find all BAMs
echo "=== ALTERNATIVE SEARCH ==="
find "$BAM_DIR" -type f -name "*.bam" | grep -E "chr1|chr19" || echo "No files found with find+grep"

# Process each BAM file - METHOD 1: Explicit listing
for chrom in chr1 chr19; do
    echo "=== PROCESSING $chrom FILES ==="
    for bamfile in "$BAM_DIR"/*"$chrom"*.bam; do
        # Skip if no files match
        [ -e "$bamfile" ] || continue
        
        # Skip temporary files
        if [[ "$bamfile" == *".tmp."* ]]; then
            continue
        fi
        
        filename=$(basename "$bamfile")
        sample=$(basename "$filename" .bam)
        
        echo "PROCESSING $filename as $chrom"
        
        # Process each position for this chromosome
        for locus in "${!LOCI[@]}"; do
            if [[ "$locus" == "$chrom:"* ]]; then
                pos=${locus#*:}
                region="${chrom}:${pos}-${pos}"
                
                # Get depth
                depth=$(samtools depth -r "$region" "$bamfile" 2>/dev/null | awk -v p="$pos" '$2==p {print $3}')
                depth=${depth:-0}
                
                echo -e "${sample}\t${chrom}\t${pos}\t${depth}" >> "$OUTPUT"
                echo " - Extracted $chrom:$pos depth: $depth"
            fi
        done
    done
done

echo "✅ Depth extraction complete. Results in $OUTPUT"