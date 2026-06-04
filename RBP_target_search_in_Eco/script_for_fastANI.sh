#!/bin/bash
set -euo pipefail

THREADS=8
BLOCK_SIZE=100

### BAKTA_DIR is the directory that you run from, while FASTANI_DIR is the location for the results, and TSV_DIR is the location for the text files of file paths to run fastANI ###
BAKTA_DIR=/home/deces123/Project/Vin_Phage/bakta
FASTANI_DIR=${BAKTA_DIR}/fastANI
TSV_DIR=${FASTANI_DIR}/tsv

### Make the result and file path directories as specified above ###
mkdir -p "${FASTANI_DIR}"
mkdir -p "${TSV_DIR}"

### Work in the BAKTA_DIR ###
cd "${BAKTA_DIR}"

### Make a full genome file path list starting from the BAKTA_DIR ###

find "$(pwd)" \
    -mindepth 2 \
    -type f \
    -name "*.fna" \
    | sort > "${FASTANI_DIR}/all_genomes.txt"

### Split the full genome file path into block files, which creates files like block_00, block_01, etc. ### 

split -l ${BLOCK_SIZE} \
      -d \
      -a 2 \
      "${FASTANI_DIR}/all_genomes.txt" \
      "${FASTANI_DIR}/block_"

blocks=(${FASTANI_DIR}/block_*)

echo "Found ${#blocks[@]} blocks"

### Run FastANI sequentially for block vs. block comparisons (block_00 vs block_00, block_00 vs block_01, etc., block_06 vs block_06 ### 

for ((i=0; i<${#blocks[@]}; i++))
do
    for ((j=i; j<${#blocks[@]}; j++))
    do

        q=${blocks[$i]}
        r=${blocks[$j]}

        qname=$(basename "$q")
        rname=$(basename "$r")

        out="${TSV_DIR}/${qname}_VS_${rname}.tsv"

        if [[ -s "$out" ]]
        then
            echo "Skipping existing $out"
            continue
        fi

        echo "Running ${qname} vs ${rname}"

        fastANI \
            --ql "$q" \
            --rl "$r" \
            -t ${THREADS} \
            -o "$out"

    done
done

### Concatenate all results ###

cat ${TSV_DIR}/block_*_VS_block_*.tsv \
    > ${TSV_DIR}/all_fastani.tsv

echo "Finished."
echo "Combined results:"
echo "${TSV_DIR}/all_fastani.tsv"