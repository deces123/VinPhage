Starting from BAKTA results, the general workflow for GWAS analyses of genes correlated with phage susceptibility in E.coli is:
- To evaluate the genetic relationship of the bacteria, I calculated the Average Nucleotide Identity (ANI) and Average Amino Acid Identity (AAI) of every sample using FastANI version 1.34 (Jain et al., 2018) and fastAAI (Gerhardt et al., 2025). 
- I used Roary version 3.13.0 (Page et al., 2015) to determine the core genome and construct a gene presence/absence matrix across all bacteria. 
- I used two different genome-wide association study (GWAS) tools, pyseer version 1.4.0 (Lees et al., 2018) and Scoary2 version 0.0.15 (Roder et al., 2024), to associate bacteria genotypes with phage host range phenotype data.
- Before the GWAS, I have to take into account population structure. I used snp-sites version 2.5.1 (Page et al., 2016) to filter out invariant sites in the core genome alignment from Roary, before using FastTree version 2.2.0 (Price et al., 2010) to build a Newick tree. 
- I used unitig-caller version 1.3.2 (Holley & Melsted, 2020) to build a population graph to extract unitigs for pyseer. I then ran pyseer with a matrix that contains binary infection values for bacteria/phage, a converted Newick tree to fit with pyseer's format, and the unitigs file.
- The gene presence/absence matrix from Roary, the Newick tree, and a matrix that contains binary infection values for bacteria/phage were used to run Scoary2.

Below, I explained what I did, generally speaking. The miscellaneous steps like setting up conda environments, ensuring files are in the correct formats, naming conventions are the same between files, etc., are up to the readers/users. This assumes you have done those steps yourself. I do not guarantee these tools will work the same if you use different versions of them. 

The RBP_target_search_based_on_keywords.py current use is mostly to confirm that certain genes/proteins are not recognized/annotated in the original BAKTA faa files (you may use it if you work with known datasets where the target of phage RBPs is already known).

- fastANI: 
The script_for_fastANI.sh is used only because a personal laptop does not have enough RAM to do an all-vs-all Average Nucleotide Identity (ANI) comparison. 
If working on a HPC with enough RAM and cores, just skip the script and go straight to command-line, where reference.txt and query.txt contains the file paths to your samples. 

(example "fastANI --ql reference.txt --rl query.txt -o fastani.out -t 10")

- fastAAI: 
Just move all the .faa files from BAKTA to the same folder of your choosing, and run the following (adjust threads used accordingly to your RAM):

"fastaai build_db --proteins fastaai_all_proteins/ --threads 6 --verbose --output fastAAI --database fastAAI_db.db --compress"
"fastaai db_query --query fastAAI/database/fastAAI_db.db --target fastAAI/database/fastAAI_db.db --threads 6 --verbose --output_style matrix --do_stdev --output fastAAI"

- Roary (Make a text file that contains paths to the gff files. If the gff files do not contain genomic sequences, add it from the FASTA file to the gff. Onwards, run like below): ("-e --mafft" to quickly generate a core gene alignment; "- p" is for thread allocation; "-i" is for minimum percentage identity for blastp)

"roary -e --mafft -i 95 -p 8 $(cat gff_cluster1_list.txt)"

- snp-sites (uses the core_gene_alignment.aln output by roary):

"snp-sites -c core_gene_alignment.aln > core_snps.aln"

- FastTree (uses the core_snps.aln produced by snp-sites):

"FastTree -fastest core_snps.aln > core_GWAS.tree"

- unitig-caller (prepare a text file that contains the file paths to all the fasta files; --kmer can be specified for the kmer size used to built the graph. By default this is 31 bp):

"unitig-caller --call --refs refs.txt --out out_prefix"

- Before pyseer, you need to convert the Newick tree produced by FastTree into a format usable by pyseer (the below flag is to make an output for use with pyseer's mixed model):
"python phylogeny_distance.py --lmm core_genome.tree > phylogeny_similarity.tsv"

- pyseer only performs GWAS per trait/phage. If you want to run pyseer automatically through multiple traits/phages, uses the script script_to_run_pyseer_recursively.sh by opening the script and edit the locations of your input files, output files, and other parameters to your discretion. Otherwise, just use command-line, and uses "pyseer --help" for instructions on how to do so.

- Scoary2 (gene_presence_absence.Rtab is returned from roary. traits.csv is a bacteria-phage infection matrix that you prepare):

"scoary2 --genes gene_presence_absence.Rtab --gene-data-type 'gene-count:\t' --traits traits.csv --outdir scoary_result --n-permut 1000"
