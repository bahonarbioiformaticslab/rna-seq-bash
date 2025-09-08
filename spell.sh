#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
THREADS=8
OUTDIR='out'

# Your SRR IDs
SRR_IDS=(
SRR9879594
SRR9879595
SRR9879596
SRR9879597
SRR9879598
SRR9879599
SRR9879600
SRR9879601
SRR9879602
SRR9879603
SRR9879604
SRR9879605
)





# === 2. Loop over SRR IDs ===
for SRR in "${SRR_IDS[@]}"; do
    echo "=== Processing $SRR ==="
    slp="${SRR}.sra"
 
    if [ ! -f "$FASTQ" ]; then
        echo "❌ FASTQ file $FASTQ not found! Skipping $SRR..."
        continue
    fi
    #do splite 
    fastq-dump --split-files $slp

    echo "✅ Finished $SRR. fastq has made"
done

# for ative file do chmode +x name_file
# for run do ./name_file

