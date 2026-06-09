#!/usr/bin/env bash

set -euo pipefail

# INPUT FILES
TRAITS="traits_fixed.csv"
KMERS="/home/deces123/Project/Vin_Phage/bakta/unitig_results.pyseer.gz"
SIM="phylogeny_distances_pyseer.tsv"

# OUTPUT DIRS
PHENO_DIR="pyseer_phenotype_input"
OUTDIR="pyseer_results"
LOGDIR="${OUTDIR}/logs"

# Phenotypes to skip
SKIP_PHENOS=()

# Get phenotype names from header (skip first column: sample)
phenotypes=$(head -1 "${TRAITS}" | tr ',' '\n' | tail -n +2)

# Loop through all phenotypes
for pheno in ${phenotypes}; do

    # Skip already completed phenotypes
    if [[ " ${SKIP_PHENOS[*]} " =~ " ${pheno} " ]]; then
        echo "[SKIP] ${pheno} already completed"
        continue
    fi

    echo "========================================"
    echo "Running pyseer for phenotype: ${pheno}"
    echo "========================================"

    # phenotype file
    PHENO_FILE="${PHENO_DIR}/${pheno}.tsv"

    # pyseer outputs
    RESULT_FILE="${OUTDIR}/${pheno}_kmers.txt"
    PATTERN_FILE="${OUTDIR}/${pheno}_kmer_patterns.txt"

    # logs
    LOG_FILE="${LOGDIR}/${pheno}.log"

    # Create phenotype TSV
    csvtk cut -f sample,"${pheno}" "${TRAITS}" \
        | tr ',' '\t' \
        > "${PHENO_FILE}"

    echo "[INFO] Created phenotype file: ${PHENO_FILE}" | tee "${LOG_FILE}"

    # Run pyseer
    pyseer \
        --lmm \
        --phenotypes "${PHENO_FILE}" \
        --kmers "${KMERS}" \
        --similarity "${SIM}" \
        --output-patterns "${PATTERN_FILE}" \
        --cpu 6 \
        > "${RESULT_FILE}" \
        2>> "${LOG_FILE}"

    echo "[INFO] Finished phenotype: ${pheno}" | tee -a "${LOG_FILE}"
    echo "" | tee -a "${LOG_FILE}"

done

echo "All pyseer runs completed."