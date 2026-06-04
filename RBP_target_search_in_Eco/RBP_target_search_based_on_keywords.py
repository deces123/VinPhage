import os
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
Expected input structure:

BASE_DIR/
|-- sample1/
|   |-- sample1.ffn
|   |-- sample1.faa
|
|-- sample2/
|   |-- sample2.ffn
|   |-- sample2.faa
"""

### Replace with your relevant paths ###
BASE_DIR = "/home/deces123/Project/Vin_Phage/bakta/"
OUTPUT_DIR = "/home/deces123/Project/Vin_Phage/rbp_target_results_bacteria"

### Keyword settings ###
RBP_STRONG = [
    "LPS", 
    "O-antigen", 
    "H-antigen", 
    "OmpC", 
    "OmpF",
    "CPS",
    "K-antigen", 
    "capsule", 
    "WTA"
    "GlcNAc",
    "Capsular Polysaccharide", 
    "LamB", 
    "TolC",
    "OprM", 
    "OmpF", 
    "OmpC",
    "Pili",
    "Flagella",
    "Flagellin"
]

RBP_WEAK = [
]

BAD_KEYWORDS = [
]

LENGTH_THRESHOLD = 500

# ------------------------------------------------------------ #

### Create output directory ###
os.makedirs(OUTPUT_DIR, exist_ok=True)


### Group all subfolders as samples ###
def group_samples(base_dir):

    samples = {}

    for folder in os.listdir(base_dir):

        folder_path = os.path.join(base_dir, folder)

        if os.path.isdir(folder_path):
            samples[folder] = folder_path

    return samples


### Find relevant file ###
def find_file(folder, extension):

    for f in os.listdir(folder):

        if f.endswith(extension):
            return os.path.join(folder, f)

    return None


### Load annotations directly from FAA FASTA headers ###
def load_bakta_annotations(faa_file):

    annotations = {}

    for record in SeqIO.parse(faa_file, "fasta"):

        gene_id = record.id

        # Full FASTA description
        desc = record.description

        # Remove duplicate ID from beginning
        product = desc.replace(gene_id, "").strip()

        annotations[gene_id] = product

    return annotations


### Load protein sequences ###
def load_proteins(faa_file):

    proteins = {}

    for record in SeqIO.parse(faa_file, "fasta"):

        proteins[record.id] = {
            "sequence": str(record.seq),
            "length": len(record.seq)
        }

    return proteins


### Assign scores ###
def score_protein(product, length):

    product_l = product.lower()

    score = 0

    # Strong RBP keywords
    if any(k.lower() in product_l for k in RBP_STRONG):
        score += 3

    # Weak keywords
    elif any(k.lower() in product_l for k in RBP_WEAK):
        score += 1

    # Length filter bonus
    if length >= LENGTH_THRESHOLD:
        score += 1

    # Hypothetical protein bonus
    if "hypothetical" in product_l:
        score += 1

    return score


### Exclude obvious non-RBPs ###
def is_bad(product):

    product_l = product.lower()

    return any(b.lower() in product_l for b in BAD_KEYWORDS)


### Extract candidate RBPs ###
def extract_rbps(sample, annotations, proteins):

    results = []

    for gene_id, product in annotations.items():

        if gene_id not in proteins:
            continue

        length = proteins[gene_id]["length"]

        score = score_protein(product, length)

        # -------- FILTERS -------- #

        if score < 2:
            continue

        if is_bad(product):
            continue

        results.append({
            "sample": sample,
            "gene_id": gene_id,
            "product": product,
            "length": length,
            "score": score
        })

    return sorted(results, key=lambda x: x["score"], reverse=True)


### Write candidate protein FASTA ###
def write_fasta(sample, results, proteins, out_path):

    records = []

    for r in results:

        gene_id = r["gene_id"]

        if gene_id not in proteins:
            continue

        seq = proteins[gene_id]["sequence"]

        rec = SeqRecord(
            Seq(seq),
            id=gene_id,
            description=f"{sample} | {r['product']} | score={r['score']}"
        )

        records.append(rec)

    if records:

        with open(out_path, "w") as out:
            SeqIO.write(records, out, "fasta")


### Main workflow ###
def main():

    samples = group_samples(BASE_DIR)

    all_results = []

    print(f"Found {len(samples)} samples\n")

    for sample, folder in samples.items():

        print(f"Processing {sample}")

        faa_file = find_file(folder, ".faa")
        ffn_file = find_file(folder, ".ffn")

        print("  faa:", faa_file)
        print("  ffn:", ffn_file)

        if faa_file is None:
            print("  Missing .faa file")
            continue

        annotations = load_bakta_annotations(faa_file)

        proteins = load_proteins(faa_file)

        print(f"  Loaded {len(annotations)} annotations")
        print(f"  Loaded {len(proteins)} proteins")

        results = extract_rbps(sample, annotations, proteins)

        print(f"  Found {len(results)} RBP candidates")

        if not results:
            continue

        ### Write FASTA ###
        fasta_out = os.path.join(
            OUTPUT_DIR,
            f"{sample}_rbp.faa"
        )

        write_fasta(sample, results, proteins, fasta_out)

        ### Write per-sample TSV ###
        out_file = os.path.join(
            OUTPUT_DIR,
            f"{sample}_rbp.tsv"
        )

        with open(out_file, "w") as out:

            writer = csv.DictWriter(
                out,
                fieldnames=results[0].keys(),
                delimiter="\t"
            )

            writer.writeheader()
            writer.writerows(results)

        all_results.extend(results)

    ### Combined output ###
    if not all_results:

        print("\nNo RBP candidates found across all samples.")
        return

    combined_file = os.path.join(
        OUTPUT_DIR,
        "all_rbp_candidates.tsv"
    )

    with open(combined_file, "w") as out:

        writer = csv.DictWriter(
            out,
            fieldnames=all_results[0].keys(),
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(all_results)

    print("\nDone!")
    print(f"Results written to: {OUTPUT_DIR}")


### Run ###
if __name__ == "__main__":
    main()