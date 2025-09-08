#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
THREADS=8
REFERENCE_INDEX="/media/bbl/201da159-b034-465f-9ef6-6e4f133a16b3/sratoolkit.3.0.7-ubuntu64/bin/MOHSENIDATA/DB
/sol_index"
GTF="GCA_000188115.5_SL4.0_genomic.gbff"
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

    FASTQ="${SRR}_1.fastq"
    SAM="$OUTDIR/${SRR}.sam"
    BAM="$OUTDIR/${SRR}.bam"
    SORTED_BAM="$OUTDIR/${SRR}.sorted.bam"
    COUNT="$OUTDIR/${SRR}.count"

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

