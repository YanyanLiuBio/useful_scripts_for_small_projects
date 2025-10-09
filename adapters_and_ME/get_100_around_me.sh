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

# Extract sequences with 100bp before and after adapter1 (FASTA format)
echo "Extracting sequences for adapter 1..."
samtools fastq "$bam_file" | \
awk -v adapter="$adapter1" -v sample="$sample_id" '
NR%4==1 { header = $0 }
NR%4==2 { 
    seq = $0
    getline; getline  # Skip quality score line and next header
    qual = $0
    pos = index(seq, adapter)
    if (pos > 0) {
        start = pos - 100
        if (start < 1) start = 1
        end = pos + length(adapter) + 99
        if (end > length(seq)) end = length(seq)
        extracted_seq = substr(seq, start, end - start + 1)
        extracted_qual = substr(qual, start, end - start + 1)
        print ">" sample "_" substr(header,2) "_adapter1_flanking"
        print extracted_seq
    }
}' | seqtk seq -A > "${sample_id}_adapter1_flanking_sequences.fasta"

# Extract sequences with 100bp before and after adapter2 (FASTA format)
echo "Extracting sequences for adapter 2..."
samtools fastq "$bam_file" | \
awk -v adapter="$adapter2" -v sample="$sample_id" '
NR%4==1 { header = $0 }
NR%4==2 { 
    seq = $0
    getline; getline  # Skip quality score line and next header
    qual = $0
    pos = index(seq, adapter)
    if (pos > 0) {
        start = pos - 100
        if (start < 1) start = 1
        end = pos + length(adapter) + 99
        if (end > length(seq)) end = length(seq)
        extracted_seq = substr(seq, start, end - start + 1)
        extracted_qual = substr(qual, start, end - start + 1)
        print ">" sample "_" substr(header,2) "_adapter2_flanking"
        print extracted_seq
    }
}' | seqtk seq -A > "${sample_id}_adapter2_flanking_sequences.fasta"

echo "Adapter 1 flanking sequences written to ${sample_id}_adapter1_flanking_sequences.fasta"
echo "Adapter 2 flanking sequences written to ${sample_id}_adapter2_flanking_sequences.fasta"
