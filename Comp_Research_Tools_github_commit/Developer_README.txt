reference_set.py is the executable 


The Snakemake workflow:
all the programs the are executed by the workflow are saved in the evo_programs subfolder. 
All the programs can executed without providing any parameters, as the default parameters are saved with argparse.

in the evo_data folder all the data that is nessicary successfully start the workflow needs to be provided. 
The nessisary files are: 

nodes.dmp
names.dmp

-from the ncbi taxonomy website: https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
-this version is from 15.12.2024

ranks_from_ncbi.tsv

-from https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/taxonomy/
-can be found under under: RankType Enumeration



species.txt
-from: https://ftp.ensemblgenomes.org/pub/release-54/species.txt

species_EnsemblVertebrates.txt
-from: https://ftp.ensembl.org/pub/release-107/species_EnsemblVertebrates.txt

Both the species.txt and species_EnsemblVertebrates.txt can be replaced by up to date versions. 
By deleting the contents of the folder evo_data/results/ 
and restarting the workflow. 

After replacing those two files the workflow should take less than 3 minutes to regenerate all the files in the evo_data/results folder, 
that are nessisary to run the reference_set.py program and create a reference_set that is fit for the specific Ensembl and Ensemblgenomes version. 

The programs of the evo_programs/ folder can also be executed manually: 

transform_ranks.py
-reads in ranks_from_ncbi.tsv , nodes.dmp
-creates a dictionary from the file: ranks_from_ncbi.tsv with numbers corresponding to the order of their appearance in the file ranks_from_ncbi.tsv
-saves this dictionary evo_data/results/ncbi_rank_dict.txt
-reads in the first 3 columns of nodes.dmp and replaces the strings from the third column "rank" with the numeric values from the ncbi_rank_dict.txt file
-outputs the first 3 columns of nodes.dmp with the numeric rank values as nodes_copy.tsv

create_ncbi_tree.py
-takes nodes_copy.tsv as input
-creates the ncbi tree structure from nodes_copy.tsv
-saves the tree to evo_data/ncbitree.pickle
-this process takes 400 000 seconds , so more than 4 days


cleanup_species.py
-reads in species.txt 
-outputs species_cleaned.txt
-this avoids errors with the UTF-8 encoding

tree_indices.py
-reads in the ncbitree.pickle
-reads in the names.dmp and saves only the names for each unique node_id to a dataframe.
-reads in the nodes_copy.tsv file
-generates the index_genomes_df which is read in by reference_set.py and used to randomly choose the node_ids of the descendants of the lowerbound set from. 







