#!/bin/bash


awk -F',' 'NR==1 {print "sample,crosstalk_ratio,largest,second_largest,crosstalk_counts,unmapped,total"; next} 
{
    # Find max1 and max2 (excluding unmapped and total columns)
    max1=0; max2=0; max1_val=0; max2_val=0
    for(i=2; i<=NF-2; i++) {
        if($i > max1_val) {
            max2_val=max1_val; 
            max1_val=$i
        }
        else if($i > max2_val) {
            max2_val=$i
        }
    }
    
    # Calculate crosstalk (sum of all except max1, max2, unmapped, total)
    crosstalk=0
    for(i=2; i<=NF-2; i++) {
        if($i != max1_val && $i != max2_val) crosstalk += $i
    }
    
    # Get unmapped and total
    unmapped=$(NF-1)
    total=$NF
    
    # Calculate ratio
    if(total > 0) {
        ratio = crosstalk/total
    } else {
        ratio = 0
    }
    
    # Accumulate for overall statistics
    total_crosstalk += crosstalk
    total_reads += total
    
    printf "%s,%.6f,%d,%d,%d,%d,%d\n", $1, ratio, max1_val, max2_val, crosstalk, unmapped, total
}
END {
    if(total_reads > 0) {
        overall_ratio = total_crosstalk/total_reads
        printf "\nOverall crosstalk ratio: %.6f\n", overall_ratio
    }
}' 20251015_ONT-minimap2_crosstalk.counts.csv