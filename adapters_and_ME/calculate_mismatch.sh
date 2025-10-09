#!/bin/bash

# Inputs
bam="LG313.TPase.Col-CC_v1.sorted.bam"
bed="centromeres.bed"

# Extract sample ID from BAM filename
sample_id=$(basename "$bam" .sorted.bam)

# Output files 
stats_csv="${sample_id}_read_stats.csv"
summary_csv="${sample_id}_read_stats_summary.csv"

# Extract reads overlapping regions in BED and compute per-read stats
samtools view -L "$bed" "$bam" \
  | awk 'BEGIN{OFS=","}
    {
      # $1=qname, $6=CIGAR, $10=SEQ, tags start at $12
      nm="NA";
      for(i=12;i<=NF;i++){
        if($i ~ /^NM:i:/){ split($i,a,":"); nm=a[3]; break }
      }
      readlen=length($10);
      if(nm=="NA") nm="";
      mismatch_rate = (nm=="" ? "" : (nm / (readlen>0 ? readlen : 1)));
      print $1, readlen, nm, mismatch_rate;
    }' > "$stats_csv"

# Compute summary statistics
awk -F, -v sample="$sample_id" 'NR>0{
  n++; 
  sumlen+=$2;
  if($3!=""){ 
    nm_sum+=($3+0); 
    nm_n++;
    total_mismatches+=($3+0);
    total_bases+=$2;
  }
}
END{
  mean_readlen = (n ? sumlen/n : 0);
  mean_NM = (nm_n ? nm_sum/nm_n : "NA");
  mismatch_pct = (total_bases > 0 ? (total_mismatches / total_bases * 100) : "NA");
  
  print "Sample=" sample, "reads=" n, "mean_readlen=" mean_readlen, "mean_NM=" mean_NM, "mismatch_pct=" mismatch_pct "%";
}' "$stats_csv" > "$summary_csv"

echo "Per-read stats written to $stats_csv"
echo "Summary stats written to $summary_csv"n


