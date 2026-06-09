#!/bin/bash

### Edit here for the file path of your input scoary folder, and your output tsv
BASE="/home/deces123/Project/Vin_Phage/bakta/scoary_Ecoli/scoary"
OUT="${BASE}/significant_genes.tsv"

header_written=0
> "$OUT"

### Edit the two 0.05 values below to your desired threshold for best_fisher_q and best_empirical_p of the scoary traits
awk -F'\t' '
NR > 1 && $3 < 0.05 && $4 < 0.05 {
    print $1
}
' "${BASE}/summary.tsv" | while read trait; do

    FILE="${BASE}/traits/${trait}/result.tsv"

    [[ -f "$FILE" ]] || continue

    if [[ $header_written -eq 0 ]]; then
        head -n 1 "$FILE" | awk -F'\t' '
        BEGIN{OFS="\t"}
        {
            print "Trait", $0
        }' > "$OUT"

        header_written=1
    fi

    awk -F'\t' -v trait="$trait" '
    BEGIN{OFS="\t"}

    NR==1 {
        for(i=1;i<=NF;i++) {
            if($i=="fisher_q") fq=i
            if($i=="empirical_p") ep=i
        }
        next
    }
### Edit the two 0.05 values below to your desired threshold for best_fisher_q and best_empirical_p of the genes in each trait
    $fq < 0.05 && $ep < 0.05 {
        print trait, $0
    }
    ' "$FILE" >> "$OUT"

done
