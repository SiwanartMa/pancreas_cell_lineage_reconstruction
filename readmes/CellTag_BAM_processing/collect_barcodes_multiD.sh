#!/bin/bash
#SBATCH --partition=msc_appbio

# Output directory for collected barcode files
OUTDIR="/scratch/users/k24023654/researchProject/ScaleRNA/outputs/celltag_barcodes"
mkdir -p "$OUTDIR"

# Process each batch (D2, D6, D8, D10)
for batch in D2 D6 D8 D10; do
    echo "Searching barcodes for $batch..."

    SOLO_BASE="/scratch/users/k24023654/researchProject/ScaleRNA/originals/alignment/CellTag-${batch}.ScaleRNA/split"

    # Loop over all solo output folders
    for dir in ${SOLO_BASE}/*.star.solo; do
        sample=$(basename "$dir" .star.solo)  # e.g., CellTag-D2.ScaleRNA.41
        barcode_file="${dir}/GeneFull_Ex50pAS/raw/barcodes.tsv"

        if [ -f "$barcode_file" ]; then
            dest="${OUTDIR}/${sample}.barcodes.tsv"
            echo "Copying $barcode_file → $dest"
            cp "$barcode_file" "$dest"
        else
            echo "Missing: $barcode_file"
        fi
    done
done

echo "All barcodes.tsv files collected in: $OUTDIR"
