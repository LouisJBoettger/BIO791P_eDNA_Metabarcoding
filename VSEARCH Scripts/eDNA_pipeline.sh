#!/bin/bash
set -euo pipefail

# Check that a barcode was supplied
if [ -z "$1" ]; then
    echo "Usage: $0 <barcode>"
    exit 1
fi

BARCODE=$1

OUT="${BARCODE}_Mammal_eDNAsh_test.txt"

# sections commented out of the run if already done by a previous pipeline and as such are not needed - reactivate as required

# Merge FASTQ files
#cat ./$BARCODE/*.fastq.gz > ${BARCODE}_merged.fastq.gz || exit 1

# Quality assessment
# NanoPlot --fastq ${BARCODE}_merged.fastq.gz --loglength || exit 1

# Quality filtering
# the line below is commented out as files were unzipped in a previous run - reactivate if needed, and edit the NanoFilt command below to cooperate
# gunzip -c ${BARCODE}_merged.fastq.gz | \
NanoFilt -q 10 --length 300 --maxlength 1500 \
    < ${BARCODE}_merged.fastq \
    > ${BARCODE}_filtered_Mammal.fastq || exit 1

# Adapter trimming
porechop \
-i ${BARCODE}_filtered_Mammal.fastq \
-o ${BARCODE}_trimmed_Mammal.fastq || exit 1

# Convert FASTQ to FASTA
vsearch \
--fastq_filter ${BARCODE}_trimmed_Mammal.fastq \
--fastaout ${BARCODE}_trimmed_Mammal.fasta \
--fastq_qmax 50 || exit 1

# Dereplicate
vsearch \
--derep_fulllength ${BARCODE}_trimmed_Mammal.fasta \
--output ${BARCODE}_dereplicated_Mammal.fasta \
--sizeout || exit 1

# Cluster OTUs
vsearch \
--cluster_size ${BARCODE}_dereplicated_Mammal.fasta \
--id 0.99 \
--centroids ${BARCODE}_otus_Mammal.fasta \
--relabel OTU_ \
--sizeout || exit 1

# Remove Chimeras
vsearch \
    --uchime_denovo ${BARCODE}_otus_Mammal.fasta \
    --nonchimeras ${BARCODE}_nochim_Mammal.fasta || exit 1

# Search against reference database
vsearch \
--usearch_global ${BARCODE}_nochim_Mammal.fasta \
--db Mammal_derep_ref.fasta \
--id 0.95 \
--blast6out "${OUT}.tmp" \
--top_hits_only \
--mincols 80 \
--maxaccepts 10 \
--maxrejects 256 || exit 1

# Add species name as final column
awk '
BEGIN {
    while ((getline < "mammal_species.txt") > 0) {
        species[$1] = $2 " " $3
    }
}
{
    acc = $2
    sub(/;size=.*/, "", acc)
    print $0 "\t" species[acc]
}
' "${OUT}.tmp" > "$OUT"

rm "${OUT}.tmp"




echo "$BARCODE completed"
