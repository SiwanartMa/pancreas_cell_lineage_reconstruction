#!/bin/bash

# Input and output directories
FASTQ_DIR="/scratch/prj/gtrm_barcoding/scratch_tmp/00_fastq"
OUT_DIR="/scratch/users/k24023654/researchProject/ScaleRNA/check_celltags"

# Output summary file
OUT_FILE="$OUT_DIR/celltag_screen_summary_detailed.tsv"

# Define regex for CellTag patterns (strict vs loose)
LOOSE_SEQ="GGT[ACGT]\\{8\\}GAATTC"
STRICT_SEQ="CCGGT[ACGT]\\{8\\}GAATTC"

# Write header
echo -e "Sample\tRead\tTotalReads\tLooseHits\tStrictHits\tLoosePct\tStrictPct" > "$OUT_FILE"

# Loop over all FASTQ files
cd "$FASTQ_DIR"
for fq in ScaleRNA-AP1-*_001.fastq.gz; do
    sample=$(basename "$fq" | cut -d'_' -f1-2)
    readtype=$(basename "$fq" | cut -d'_' -f2)  # R1 or R2

    # Count total number of reads (FASTQ = 4 lines per read)
    total_reads=$(zcat "$fq" | wc -l)
    total_reads=$((total_reads / 4))

    # Count CellTag hits
    loose_hits=$(zcat "$fq" | grep -o "$LOOSE_SEQ" | wc -l)
    strict_hits=$(zcat "$fq" | grep -o "$STRICT_SEQ" | wc -l)

    # Calculate percentages
    if [ "$total_reads" -gt 0 ]; then
        loose_pct=$(awk "BEGIN {printf \"%.4f\", ($loose_hits / $total_reads) * 100}")
        strict_pct=$(awk "BEGIN {printf \"%.4f\", ($strict_hits / $total_reads) * 100}")
    else
        loose_pct="NA"
        strict_pct="NA"
    fi

    # Output results to file
    echo -e "${sample}\t${readtype}\t${total_reads}\t${loose_hits}\t${strict_hits}\t${loose_pct}\t${strict_pct}" >> "$OUT_FILE"
done

echo "Finished scanning CellTag motifs. Summary saved to:"
echo "$OUT_FILE"
