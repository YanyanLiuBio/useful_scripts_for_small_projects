## WGS Fold Penalty Calculator - Multi-file Version

This python script process multiple Picard files at once and automatically
matches them to samples in the mosdepth file.

Usage:
    python calculate_wgs_fold_penalties_multi.py \
        -m mosdepth.tsv \
        -p file1.wgs.txt file2.wgs.txt file3.wgs.txt \
        --excel results.xlsx