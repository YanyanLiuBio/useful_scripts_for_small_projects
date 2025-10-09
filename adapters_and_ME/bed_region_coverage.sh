#!/bin/bash
set -euo pipefail

# Usage: ./arms_coverage.sh <reference.fa> <bam> <centromeres.bed>

ref_fa="$1"
bam="$2"
centromeres_bed="$3"

# extract sample ID from BAM filename
sample_id=$(basename "$bam" .bam)

# create chrom.sizes from reference FASTA index
cut -f1,2 "${ref_fa}.fai" > chrom.sizes

# create arms bed (whole chrom minus centromeres)
bedtools subtract -a <(awk '{print $1"\t0\t"$2}' chrom.sizes) -b "$centromeres_bed" > arms.bed

# compute mean coverage per arms interval
bedtools coverage -a arms.bed -b "$bam" -mean > "${sample_id}.arms.coverage.mean.tsv"

# compute overall mean across arms
awk -v sample="$sample_id" '{sum+=$7; n++} END{printf "%s,%d,%.2f\n", sample, n, sum/n}' "${sample_id}.arms.coverage.mean.tsv" \
    > "${sample_id}.arms.coverage.summary.csv"

echo "Done! Output files:"
echo "  Interval-level coverage: ${sample_id}.arms.coverage.mean.tsv"
echo "  Overall mean coverage: ${sample_id}.arms.coverage.summary.csv"
