#!/bin/bash

# Define input/output directories and scripts
BASE_INPUT_DIR="/scratch/users/k24023654/researchProject/ScaleRNA/originals/alignment"
SCRIPTS_DIR="/scratch/users/k24023654/researchProject/CellTag_practice/script"
OUTPUT_DIR="./results_v1_motif"
V1_MOTIF="GGT[ACTG]{8}GAATTC"  # v1 CellTag motif
REGEX="CCGGT([ACTG]{8})GAATTC"  # stricter v1 CellTag motif
SAMPLES=("CellTag-D2.ScaleRNA" "CellTag-D6.ScaleRNA" "CellTag-D8.ScaleRNA" "CellTag-D10.ScaleRNA" "CellTag-iPSCs-Frozen.ScaleRNA")

# Load tools
module load samtools/1.17-gcc-13.2.0-python-3.11.6
module load r/4.3.0-gcc-13.2.0-withx-rmath-standalone-python-3.11.6

# Create output directory
mkdir -p "$OUTPUT_DIR"

# -----------------------------
# Main loop over all CellTag sample batches
# -----------------------------
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_DIR="${BASE_INPUT_DIR}/${SAMPLE}/split"

    for BAM in "$SAMPLE_DIR"/*.star.align/*.bam; do
        ALIGN_DIR=$(basename "$(dirname "$BAM")")
        SAMPLE_ID=${ALIGN_DIR/.star.align/} 

        echo "Processing $SAMPLE_ID"

        # Create output directory for this sample
        SAMPLE_OUT_DIR="${OUTPUT_DIR}/${SAMPLE_ID}"
        mkdir -p "$SAMPLE_OUT_DIR"
        cd "$SAMPLE_OUT_DIR"

        # Step 1: Extract reads with CellTags
        READS_OUT="${SAMPLE_ID}.celltag.reads.out"
        echo "Extracting CellTag-containing reads"
        samtools view "$BAM" | grep -P "$V1_MOTIF" > "$READS_OUT"

        # Step 2: Parse extracted reads
        PARSED_OUT="${SAMPLE_ID}.celltag.parsed.tsv"
        echo "Parsing the reads"
        "${SCRIPTS_DIR}/celltag.parse.reads.10x.sh" -v tagregex="$REGEX" "$READS_OUT" > "$PARSED_OUT"

        # Step 3: Quantify CellTags using barcodes.tsv from STARsolo
        BARCODES_FILE="${SAMPLE_DIR}/${SAMPLE_ID}.star.solo/GeneFull_Ex50pAS/raw/barcodes.tsv"
        echo "Quantifying CellTags per barcode"
        Rscript "${SCRIPTS_DIR}/matrix.count.celltags.R" "$BARCODES_FILE" "$PARSED_OUT" "$SAMPLE_ID"

        # Step 4: Return to original directory
        cd - > /dev/null
    done
done

