#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
THREADS=8
REFERENCE_INDEX="/media/bbl/201da159-b034-465f-9ef6-6e4f133a16b3/soil_s/dbtomato/solh"
GTF="Solanum_lycopersicum.SL3.0.54.gtf"
OUTDIR='out'

# Your SRR IDs
SRR_IDS=(
SRR22245988
SRR22245989
SRR22245990
SRR22245991
SRR22245992
SRR22245993

)





# === 2. Loop over SRR IDs ===
for SRR in "${SRR_IDS[@]}"; do
    echo "=== Processing $SRR ==="

    FASTQ1="${SRR}_1.fastq"
    FASTQ2="${SRR}_2.fastq"
    SAM="$OUTDIR/${SRR}.sam"
    BAM="$OUTDIR/${SRR}.bam"
    SORTED_BAM="$OUTDIR/${SRR}.sorted.bam"
    COUNT="$OUTDIR/count.${SRR}.txt"

    if [ ! -f "$FASTQ1" ]; then
        echo "❌ FASTQ file $FASTQ1 not found! Skipping $SRR..."
        continue
    fi

    # Run HISAT2
    hisat2 -p "$THREADS"  -x "$REFERENCE_INDEX"  -1 "$FASTQ1" -2 "$FASTQ2" -S "$SAM"

    # Convert SAM → BAM
    samtools view -b "$SAM" > "$BAM"

    # Sort BAM
    samtools sort "$BAM" -o "$SORTED_BAM"

    # Count reads
    htseq-count -f bam -s reverse "$SORTED_BAM" "$GTF" > "$COUNT"

    echo "✅ Finished $SRR. Counts saved to $COUNT"
done

