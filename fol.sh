#!/usr/bin/env bash

# Define the list of SRR IDs
srr_ids=(
    SRR14118679
    SRR14118680
    SRR14118681
    SRR14118682
    SRR14118683
    SRR14118684
    SRR14118685
    SRR14118686
    SRR14118687
    SRR14118688
    SRR14118689
    SRR14118690
    SRR14118691
    SRR14118692
    SRR14118693
    SRR14118694
    SRR14118695
    SRR14118696
)

# Create output directory if it doesn't exist
mkdir -p data/raw_sra
mkdir -p logs

# Loop over the list and run prefetch
for i in "${srr_ids[@]}"; do
    echo "Downloading $i..."
    prefetch "$i" -O data/raw_sra/ > logs/${i}.log 2>&1
done

