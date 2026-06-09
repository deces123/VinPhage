import pandas as pd
from pathlib import Path

### Alter this to the file path of the significant_genes.tsv from scoary
scoary = pd.read_csv(
    "/home/deces123/Project/Vin_Phage/bakta/scoary_Ecoli/scoary/significant_genes.tsv",
    sep="\t"
)

### Alter this to the file path of the step #11 results
gene_hits_dir = Path(
    "/home/deces123/Project/Vin_Phage/bakta/scoary_Ecoli/pyseer_results/bonferroni_filtered_results/gene_hits_summary"
)

matches = []

for _, row in scoary.iterrows():

    trait = str(row["Trait"]).strip()
    gene = str(row["Gene"]).strip()

    ### Alter this to the general names/extensions of the result files of step #11
    pyseer_file = gene_hits_dir / f"{trait}_gene_summary.txt"

    if not pyseer_file.exists():
        continue

    pyseer = pd.read_csv(pyseer_file, sep="\t")

    # Normalize gene names
    pyseer["gene"] = pyseer["gene"].astype(str).str.strip()

    hit = pyseer[pyseer["gene"] == gene]

    if hit.empty:
        continue

    for _, h in hit.iterrows():

        # Copy ALL Scoary columns
        result = row.to_dict()

        # Add Pyseer statistics
        result.update({
            "pyseer_hits": h["hits"],
            "pyseer_maxp": h["maxp"],
            "pyseer_avg_af": h["avg_af"],
            "pyseer_avg_maf": h["avg_maf"],
            "pyseer_avg_beta": h["avg_beta"]
        })

        matches.append(result)

matches_df = pd.DataFrame(matches)

### Edit the below for the file path of the result
matches_df.to_csv(
    "/home/deces123/Project/Vin_Phage/bakta/scoary_Ecoli/scoary_pyseer_overlaps.tsv",
    sep="\t",
    index=False
)

print(f"Found {len(matches_df)} matching gene-trait pairs")