#!/bin/bash

bam_file="LG313.TPase.Col-CC_v1.sorted.bam"

# Extract sample ID from BAM filename
sample_id=$(basename "$bam_file" .bam)

# Use sample ID in report filename
report_file="${sample_id}_adapter_report.txt"

adapter1="AGATGTGTATAAGAGACAG"
adapter2="CTGTCTCTTATACACATCT"

# Count reads containing each adapter
count1=$(samtools fastq "$bam_file" | awk 'NR%4==2' | grep -F -c "$adapter1")
count2=$(samtools fastq "$bam_file" | awk 'NR%4==2' | grep -F -c "$adapter2")

# Write report with sample ID
{
  echo "Adapter sequence report for sample: $sample_id"
  echo "BAM file: $bam_file"
  echo "--------------------------------------"
  echo "Adapter 1 ($adapter1): $count1 reads"
  echo "Adapter 2 ($adapter2): $count2 reads"
  echo "Total reads with adapters: $((count1 + count2))"
} > "$report_file"

echo "Report written to $report_file"
