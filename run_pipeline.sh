#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
THREADS=8
REFERENCE_INDEX="/media/bbl/201da159-b034-465f-9ef6-6e4f133a16b3/soil_s/dbtomato/solh"
GTF="Solanum_lycopersicum.SL3.0.54.gtf"
OUTDIR='out'

# Your SRR IDs
SRR_IDS=(
    SRR10550643
    SRR10550646
    SRR10550647
    SRR10550648
    SRR10550654
    SRR10550655
)





# === 2. Loop over SRR IDs ===
for SRR in "${SRR_IDS[@]}"; do
    echo "=== Processing $SRR ==="

    FASTQ="${SRR}_1.fastq"
    SAM="$OUTDIR/${SRR}.sam"
    BAM="$OUTDIR/${SRR}.bam"
    SORTED_BAM="$OUTDIR/${SRR}.sorted.bam"
    COUNT="$OUTDIR/count.${SRR}.txt"

    if [ ! -f "$FASTQ" ]; then
        echo "❌ FASTQ file $FASTQ not found! Skipping $SRR..."
        continue
    fi

    # Run HISAT2
    hisat2 -p "$THREADS" -x "$REFERENCE_INDEX" -U "$FASTQ" -S "$SAM"

    # Convert SAM → BAM
    samtools view -b "$SAM" > "$BAM"

    # Sort BAM
    samtools sort "$BAM" -o "$SORTED_BAM"

    # Count reads
    htseq-count -f bam -s reverse "$SORTED_BAM" "$GTF" > "$COUNT"

    echo "✅ Finished $SRR. Counts saved to $COUNT"
done

