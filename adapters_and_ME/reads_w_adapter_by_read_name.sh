# Get total read count for reference
total_reads=$(samtools view -c LG313.SPK3.Col-CC_v1.sorted.bam)
echo "Total reads: $total_reads"

# Get read IDs for each sequence
samtools view LG313.SPK3.Col-CC_v1.sorted.bam | awk '/CGCCGTATCATT/ {print $1}' | sort -u > adapter1_reads.txt
samtools view LG313.SPK3.Col-CC_v1.sorted.bam | awk '/AATGATACGGCG/ {print $1}' | sort -u > adapter1_RC_reads.txt

# Count unique reads for each sequence
adapter1_count=$(wc -l < adapter1_reads.txt)
adapter1_RC_count=$(wc -l < adapter1_RC_reads.txt)

echo "Unique reads with adapter1: $adapter1_count"
echo "Unique reads with adapter1_RC: $adapter1_RC_count"

# Find reads that have BOTH sequences (intersection)
comm -12 <(sort adapter1_reads.txt) <(sort adapter1_RC_reads.txt) > both_sequences_reads.txt
both_count=$(wc -l < both_sequences_reads.txt)
echo "Unique reads with BOTH sequences: $both_count"

# Find reads that have ONLY adapter1
comm -23 <(sort adapter1_reads.txt) <(sort adapter1_RC_reads.txt) > only_adapter1_reads.txt
only_adapter1_count=$(wc -l < only_adapter1_reads.txt)
echo "Unique reads with ONLY adapter1: $only_adapter1_count"

# Find reads that have ONLY adapter1_RC
comm -13 <(sort adapter1_reads.txt) <(sort adapter1_RC_reads.txt) > only_adapter1_RC_reads.txt
only_adapter1_RC_count=$(wc -l < only_adapter1_RC_reads.txt)
echo "Unique reads with ONLY adapter1_RC: $only_adapter1_RC_count"

# Summary
echo "=== SUMMARY ==="
echo "Total reads in BAM: $total_reads"
echo "Reads with adapter1 (unique): $adapter1_count"
echo "Reads with adapter1_RC (unique): $adapter1_RC_count"
echo "Reads with BOTH sequences: $both_count"
echo "Reads with ONLY adapter1: $only_adapter1_count"
echo "Reads with ONLY adapter1_RC: $only_adapter1_RC_count"
