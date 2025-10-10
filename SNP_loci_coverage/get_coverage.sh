#!/bin/bash

BAM_DIR="bam_subset"
OUT_DIR="downsampled_bams"
mkdir -p "$OUT_DIR"
TARGET_COVERAGES=(1 2 5 10 15)

# Function to format fraction correctly for samtools
format_fraction() {
    printf "%.6f" "$1" | sed 's/0*$//' | sed 's/\.$//'
}

for sample in 24300FL-07-01-01_S1_L006 24300FL-07-01-02_S2_L006; do
    for chr in 1 19; do
        input_bam="$BAM_DIR/${sample}.md_chr${chr}.bam"
        
        if [[ ! -f "$input_bam" ]]; then
            echo "File $input_bam not found, skipping..."
            continue
        fi
        
        echo "Processing $input_bam..."
        current_coverage=$(samtools depth "$input_bam" | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
        echo "Current coverage: $current_coverage"
        
        for target_cov in "${TARGET_COVERAGES[@]}"; do
            output_bam="$OUT_DIR/${sample}.md_chr${chr}.${target_cov}x.bam"
            fraction=$(awk -v c=$current_coverage -v t=$target_cov 'BEGIN {print (t>=c) ? 1.0 : t/c}')
            formatted_frac=$(format_fraction $fraction)
            
            if (( $(echo "$fraction >= 1.0" | bc -l) )); then
                echo "Creating symlink to original for ${target_cov}x coverage"
                ln -sf "../$input_bam" "$output_bam"
                ln -sf "../$input_bam.bai" "$output_bam.bai"
            else
                echo "Downsampling to ${target_cov}x (fraction: $formatted_frac)"
                seed=123  # Fixed seed for reproducibility
                samtools view -b -s $seed$formatted_frac "$input_bam" | samtools sort -o "$output_bam"
                samtools index "$output_bam"
            fi
            
            # Verify coverage
            final_cov=$(samtools depth "$output_bam" | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
            echo "Final coverage: $final_cov"
        done
    done
done