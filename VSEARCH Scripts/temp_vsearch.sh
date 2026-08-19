#!/bin/bash
set -euo pipefail

# Check that a barcode was supplied
if [ -z "$1" ]; then
    echo "Usage: $0 <barcode>"
    exit 1
fi

BARCODE=$1

OUT="${BARCODE}_All_eDNAsh_test.txt"

# Search against reference database
vsearch \
--usearch_global ${BARCODE}_nochim_Mammal.fasta \
--db COI_ref_derep_species.fasta \
--id 0.95 \
--blast6out "${OUT}.tmp" \
--top_hits_only \
--mincols 80 \
--maxaccepts 10 \
--maxrejects 256 || exit 1

# Add species name as final column
awk '
BEGIN {
    while ((getline < "accession_species_lookup.txt") > 0) {
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

