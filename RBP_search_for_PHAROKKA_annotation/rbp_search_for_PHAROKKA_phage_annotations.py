import os
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
Expected input structure:

BASE_DIR/
|-- sample1-pharokka/
|   |-- sample1.ffn
|   |-- sample1.faa
|
|-- sample2-pharokka/
|   |-- sample2.ffn
|   |-- sample2.faa
"""

### Replace with the relevant parameters, such as destination folder, which keyword to prioritize, etc. ###
BASE_DIR = "/home/deces123/Project/Vin_Phage/PROCESSING/genome-annotate/"                     
OUTPUT_DIR = "/home/deces123/Project/Vin_Phage/rbp_results"

RBP_STRONG = ["tail fiber", "tail spike", "receptor", "adhesion"]
RBP_WEAK = ["tail"]

BAD_KEYWORDS = [
    "polymerase", "ribosomal", "capsid",
    "portal", "terminase", "helicase"
]

LENGTH_THRESHOLD = 500
# -------------------------- #

os.makedirs(OUTPUT_DIR, exist_ok=True) ### Make new folder, or if it exists, then ignores.


### Group samples; only works if the folder structure consists of multiple subfolders with the word "pharokka" ###
def group_samples(base_dir):
    samples = {}

    for folder in os.listdir(base_dir):
        if folder.endswith("-pharokka"):
            sample = folder.replace("-pharokka", "")
            samples[sample] = os.path.join(base_dir, folder)

    return samples


### Get the relevant file paths ###
def find_file(folder, extension):
    for f in os.listdir(folder):
        if f.endswith(extension):
            return os.path.join(folder, f)
    return None


### Load Pharokka functional annotations from a .functions file. ###
def load_pharokka_functions(func_file):
    annotations = {}

    with open(func_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue

            sample, gene_id, product = parts[:3]
            annotations[gene_id] = product

    return annotations


### Load protein sequences from a FASTA (.faa) files ###
def load_proteins(faa_file):
    proteins = {}

    for record in SeqIO.parse(faa_file, "fasta"):
        proteins[record.id] = {
            "sequence": str(record.seq),
            "length": len(record.seq)
        }

    return proteins


### Assigning an arbitrary score for each grepped row (only exists since I was unsure if any protein that contain the word "tail" automatically has to do with "RBPs", also because unsure if I was supposed to find a way to further annotate the hypothetical proteins) ###
def score_protein(product, length):
    product_l = product.lower()
    score = 0

    # strong RBP keywords
    if any(k in product_l for k in RBP_STRONG):
        score += 3

    # weak (generic tail)
    elif any(k in product_l for k in RBP_WEAK):
        score += 1

    # length (exists because BLAST may potentially give a high identity for very short sequence)
    if length >= LENGTH_THRESHOLD:
        score += 1

    # hypothetical bonus
    if "hypothetical" in product_l:
        score += 1

    return score


### Exclude irrelevant annotations ###
def is_bad(product):
    product_l = product.lower()
    return any(b in product_l for b in BAD_KEYWORDS)


### Identify candidate RBPs from the score_protein function above and return the results ###
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


### Running the above functions while printing out certain warnings if certain scenarios were to occur ###
def main():
    samples = group_samples(BASE_DIR)
    all_results = []

    print(f"Found {len(samples)} samples\n")

    for sample, folder in samples.items():
        print(f"Processing {sample}")

        func_file = find_file(folder, ".functions")
        faa_file = find_file(folder, ".faa")

        print("  functions:", func_file)
        print("  faa:", faa_file)

        if func_file is None or faa_file is None:
            print("  Missing required files")
            continue

        annotations = load_pharokka_functions(func_file)
        proteins = load_proteins(faa_file)

        print(f"  Loaded {len(annotations)} annotations")
        print(f"  Loaded {len(proteins)} proteins")

        results = extract_rbps(sample, annotations, proteins)
        fasta_out = os.path.join(OUTPUT_DIR, f"{sample}_rbp.faa")
        write_fasta(sample, results, proteins, fasta_out)

        print(f"  Found {len(results)} RBP candidates")

        if not results:
            continue

        # Write per-sample candidate table (.tsv) #
        out_file = os.path.join(OUTPUT_DIR, f"{sample}_rbp.tsv")
        with open(out_file, "w") as out:
            writer = csv.DictWriter(out, fieldnames=results[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(results)

        all_results.extend(results)

    # combined output
    if not all_results:
        print("\nNo RBP candidates found across all samples.")
        return

    combined_file = os.path.join(OUTPUT_DIR, "all_rbp_candidates.tsv")
    with open(combined_file, "w") as out:
        writer = csv.DictWriter(out, fieldnames=all_results[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(all_results)

    print("\nDone!")
    print(f"Results written to: {OUTPUT_DIR}")

### Write candidate RBP protein sequences to FASTA format. ###
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

### Just run with "python test.py"/Remember to load the environment that has python ###
if __name__ == "__main__":
    main()
