#!/bin/bash

# Directory containing BAM files
BAM_DIR="bam"

# Output directory
OUT_DIR="extracted_chromosomes"
mkdir -p "$OUT_DIR"

# List of BAM files to process
BAM_FILES=(
    "24300FL-07-01-01_S1_L006.md.bam"
    "24300FL-07-01-02_S2_L006.md.bam"
)

# Chromosomes to extract
CHROMOSOMES=("chr1" "chr19")

# Process each BAM file
for bam_file in "${BAM_FILES[@]}"; do
    # Check if file exists
    if [[ ! -f "$BAM_DIR/$bam_file" ]]; then
        echo "File $BAM_DIR/$bam_file not found, skipping..."
        continue
    fi
    
    # Get base name without extension
    base_name=$(basename "$bam_file" .bam)
    
    # Process each chromosome
    for chr in "${CHROMOSOMES[@]}"; do
        output_file="$OUT_DIR/${base_name}_${chr}.bam"
        
        echo "Extracting $chr from $bam_file..."
        
        # Use samtools to extract the chromosome
        samtools view -b "$BAM_DIR/$bam_file" "$chr" > "$output_file"
        
        # Index the output BAM file
        samtools index "$output_file"
        
        echo "Created $output_file and its index"
    done
done

echo "Extraction complete!"