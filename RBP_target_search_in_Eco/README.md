Starting from BAKTA results, the general workflow is:
- To evaluate the genetic relationship of the bacteria, we calculated the Average Nucleotide Identity (ANI) and Average Amino Acid Identity (AAI) of every sample using FastANI version 1.34 (Jain et al., 2018) and fastAAI (Gerhardt et al., 2025). 
- We used Roary version 3.13.0 (Page et al., 2015) to determine the core genome and construct a gene presence/absence matrix across all bacteria. 
- We used two different genome-wide association study (GWAS) tools, pyseer version 1.4.0 (Lees et al., 2018) and Scoary2 version 0.0.15 (Roder et al., 2024), to associate bacteria genotypes with phage host range phenotype data. 
- We used unitig-caller version 1.3.2 (Holley & Melsted, 2020) to build a population graph to extract unitigs for pyseer.
- The gene presence/absence matrix from Roary was used to run scoary2.

The RBP_target_search_based_on_keywords.py current use is mostly to confirm that certain genes/proteins are not recognized/annotated in the original BAKTA faa files (you may use it if you work with known datasets where the target of phage RBPs is already known).

fastANI: 
The script_for_fastANI.sh is used only because a personal laptop does not have enough RAM to do an all-vs-all Average Nucleotide Identity (ANI) comparison. 
If working on a HPC with enough RAM and cores, just skip the script and go straight to command line (example fastANI --ql reference.txt --rl query.txt -o fastani.out -t 10)

fastAAI: 
Just move all the .faa files from BAKTA to the same folder of your choosing, and run the following (adjust threads used accordingly to your RAM):
fastaai build_db --proteins fastaai_all_proteins/ --threads 6 --verbose --output fastAAI --database fastAAI_db.db --compress
fastaai db_query --query fastAAI/database/fastAAI_db.db --target fastAAI/database/fastAAI_db.db --threads 6 --verbose --output_style matrix --do_stdev --output fastAAI

roary (Make a text file that contains paths to the gff files. If the gff files do not contain genomic sequences, add it from the FASTA file to the gff. Onwards, run like below): ("-e --mafft" to quickly generate a core gene alignment; "- p" is for thread allocation; "-i" is for minimum percentage identity for blastp)
roary -e --mafft -i 95 -p 8 $(cat gff_cluster1_list.txt)

scoary2 (gene_presence_absence.Rtab is returned from roary. traits.csv is a bacteria-phage infection matrix that you prepare):
scoary2 --genes gene_presence_absence.Rtab --gene-data-type 'gene-count:\t' --traits traits.csv --outdir scoary_result --n-permut 1000
