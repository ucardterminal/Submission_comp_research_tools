the workflow manages the following programs :


transform_ranks.py

reads in:
the nodes.dmp and the ranks_from_ncbi.tsv file from the folder evo_data/

it create pandas table of the nodes in nodes.dmp and saves them with the ranks transformed to numerical values.

it saves this table at evo_data/results/nodes_copy
the dictionary to transform these rank back gets set to lowercase with spaces inbetween and saved at:
evo_data/results/ncbi_rank_dict.tsv



taxonomy_smalltree.py

reads in:
the pandas table with the nodes and the transformed ranks from evo_data/results/nodes_copy

creates a small version of this tree to be able to create it in the scope of this project.

saves that tree to evo_data/results/ncbitree.pickle



cleanup_species.py 

reads in:
the species.txt file from evo_data/species.txt

corrects all the utf8 errors.
saves it at evo_data/results/species_cleaned.txt



tree_indices.py 
reads in:
the created tree ncbitree.pickle from : evo_data/results/ncbitree.pickle

the created evo_data/results/species_cleaned.txt file and aonther file: evo_data/species_EnsemblVertebrates.txt


the created nodes pandas table :evo_data/results/nodes_copy.tsv

and the names file for these nodes: evo_data/names.dmp

saves for each node which nodes with which names and node id lie under this node and what their names are and saves it to the tabe nodes_sum

saves the table to nodes_sum.tsv


