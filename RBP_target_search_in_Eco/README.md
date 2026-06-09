Starting from BAKTA results, the general workflow for GWAS analyses of genes correlated with phage susceptibility in E.coli is:
- To evaluate the genetic relationship of the bacteria, I calculated the Average Nucleotide Identity (ANI) and Average Amino Acid Identity (AAI) of every sample using FastANI version 1.34 (Jain et al., 2018) and fastAAI (Gerhardt et al., 2025). 
- I used Roary version 3.13.0 (Page et al., 2015) to determine the core genome and construct a gene presence/absence matrix across all bacteria. 
- I used two different genome-wide association study (GWAS) tools, pyseer version 1.4.0 (Lees et al., 2018) and Scoary2 version 0.0.15 (Roder et al., 2024), to associate bacteria genotypes with phage host range phenotype data.
- Before the GWAS, I have to take into account population structure. I used snp-sites version 2.5.1 (Page et al., 2016) to filter out invariant sites in the core genome alignment from Roary, before using FastTree version 2.2.0 (Price et al., 2010) to build a Newick tree. 
- I used unitig-caller version 1.3.2 (Holley & Melsted, 2020) to build a population graph to extract unitigs for pyseer. I then ran pyseer with a matrix that contains binary infection values for bacteria/phage, a converted Newick tree to fit with pyseer's format, and the unitigs file.
- The gene presence/absence matrix from Roary, the Newick tree, and a matrix that contains binary infection values for bacteria/phage were used to run Scoary2.

I explained what I did, generally speaking, in the Workflow_description.txt, since a lot of what I did is command lines only, so there is very little automation. The miscellaneous steps like setting up conda environments, ensuring files are in the correct formats, naming conventions are the same between files, etc., are up to the readers/users. This assumes you have done those steps yourself. I do not guarantee these tools will work the same if you use different versions of them. 

DISCLAIMER: I do not make the following scripts (they are found at https://github.com/mgalardini/pyseer/tree/master/scripts), and inside contain information on who makes them:
- count_patterns.py
- phylogeny_distance.py
- qq_plot.py
- summarise_annotations.py
