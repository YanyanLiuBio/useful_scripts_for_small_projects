#!/bin/bash

bam="LG313.TPase.Col-CC_v1.sorted.bam"
out="barcode_read_length_stats.csv"

# extract sample ID from BAM filename
sample_id=$(basename "$bam" .bam)

# Write header
echo "Sample,Barcode,Count,Min,Max,Mean,Median" > $out

for barcode in AATGATACGGCG CGCCGTATCATT; do
    txt="${barcode}_reads.txt"
    
    # Debug: check what's producing length 1
    echo "Checking short reads for $barcode..."
    samtools view "$bam" \
      | grep -F -f "$txt" \
      | awk '{len=length($10); if(len<=10) print "Read:",$1,"Seq:",$10,"Len:",len}' \
      | head -5
    
    # Main processing with validation
    samtools view "$bam" \
      | grep -F -f "$txt" \
      | awk '{
          seq=$10
          len=length(seq)
          if(len>1 && seq!="*" && seq~/^[ACGTN]+$/) print len
        }' \
      | sort -n \
      | awk -v bc="$barcode" -v sid="$sample_id" '{
          a[NR]=$1; sum+=$1
        }
        END {
          if (NR>0) {
            q2 = (NR%2 ? a[(NR+1)/2] : (a[NR/2]+a[NR/2+1])/2)
            printf("%s,%s,%d,%d,%d,%.2f,%.2f\n",
                   sid, bc, NR, a[1], a[NR], sum/NR, q2)
          }
        }' >> $out
done
