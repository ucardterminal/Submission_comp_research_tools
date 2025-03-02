# -*- coding: utf-8 -*-
import pandas as pd
import pickle
import time
import numpy as np
import csv
import argparse
from bigtree import Node, print_tree
from bigtree import Node, tree_to_dot, postorder_iter
from bigtree import Node, findall, find_children ,add_dataframe_to_tree_by_path
from bigtree import Node, find, find_name, find_path, find_relative_path, find_full_path, find_attr

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='this program creates a table with the sum of genomes for each node in the tree')
    
    # Add arguments with default values
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", default='../evo_data/ncbitree.pickle',help="filename of input file")
    parser.add_argument("--EnsemblVertebrates", default='../evo_data/species_EnsemblVertebrates.txt',help="filename of input file")
    parser.add_argument("--species_cleaned", default='../evo_data/results/species_cleaned.txt',help="filename of input file")
    parser.add_argument("--nodes", default='../evo_data/results/nodes_copy.tsv',help="filename of input file")
    parser.add_argument("--names", default='../evo_data/names.dmp',help="filename of input file")
    parser.add_argument("--input4", default='../evo_data/species_EnsemblVertebrates.txt',help="filename of input file")
    parser.add_argument("--index_genomes",default='../evo_data/results/index_genomes_df.tsv',  help="filename of nodes output file")
    parser.add_argument("--summed_nodes",default='../evo_data/results/nodes_sum_table.tsv',  help="filename of dict output file")
    args = parser.parse_args()



    with open(args.tree, 'rb') as f:  #takes 3 seconds to load ncbitree.pickle
        loadedtree = pickle.load(f)


    genomes_df = pd.read_table(args.EnsemblVertebrates, usecols=[3,5]) #'species_EnsemblVertebrates.txt'
    genomes_df.columns = ['node_id', 'assembly_accession']
    print(genomes_df)   

    species_genomes_df = pd.read_table(args.species_cleaned, usecols=[3,5],encoding='utf-8') #'species_cleaned.txt'
    species_genomes_df.columns = ['node_id', 'assembly_accession']
    print(species_genomes_df)  

    allnodes_df = pd.concat([species_genomes_df, genomes_df], ignore_index=True) # add both dataframes to each other
    print(allnodes_df) 
    print("all_nodes_df_indexed")
    print(allnodes_df)
    allnodes_df['node_id'] = allnodes_df['node_id'].astype(str)


    genomes_df = allnodes_df['node_id'].value_counts().reset_index() #create dataframe with column 1 node_id with genomes and column2 number of genomes for that node_id
    genomes_df.columns = ['node_id', 'at_node'] 
    #genomes_df.to_csv('../evo_data/genomes_df.tsv', sep='\t', index=False) #'genomes_df.tsv'

    print("genomes_df")
    print(genomes_df['at_node'].max())
    print(genomes_df)
    genomes_df['node_id'] = genomes_df['node_id'].astype(str)
    #print type of at_node
    print("type of node_id genomes_df")
    print((genomes_df['node_id'].dtype))
    print(type(genomes_df.iloc[0, 0]))

    nodes_df = pd.read_csv(args.nodes, delimiter='\t',header=None) #read in nodes table from nodes_copy.tsv
    nodes_df.columns = ['node_id', 'parent_id', 'rank']
    nodes_df['at_node'] = 0

    #print type of at_node
    nodes_df['node_id'] = nodes_df['node_id'].astype(str)
    print("type of node_id nodes_df")
    print((nodes_df['node_id'].dtype))
    print(nodes_df.iloc[0, 0])
    print("nodes_df")
    print(nodes_df)


    index_genomes_df = allnodes_df[allnodes_df['node_id'].isin(nodes_df['node_id'])] #filter out the rows in allnodes_df that are not in nodes_df
    length = len(index_genomes_df)
    print("length index_df",length)
    print("index_genomes_df")
    print(index_genomes_df)
    index_genomes_df = index_genomes_df.copy()
    index_genomes_df.loc[:, 'names'] = None

    


    data = []
    with open(args.names, 'r') as file:  #open 'names.dmp' and put it into a list to then convert it to a pandas dataframe
        reader = csv.reader(file, delimiter='|', skipinitialspace=True)
        for row in reader:
            # Strip whitespace from each cell
            row = [cell.strip() for cell in row]
            data.append(row)
        #print(data)

    # Convert to DataFrame
    names_df = pd.DataFrame(data, columns=['node_id', 'names', 'additional_info', 'specificity',"somethingelse"])

    # Drop rows where all values are NaN (if any)
    names_df.dropna(how='all', inplace=True)

    # Display the DataFrame
    print(names_df)


    print("names_df just read in ")
    print(names_df)
    names_df = names_df[names_df["specificity"] == "scientific name"] #only keep rows where the 'extra' column is 's
    print("names_df just scientific name")
    print(names_df)
    names_df = names_df[["node_id", "names"]] #keep only the node_id and scientific_name columns
    print("names_df only node_id and names")
    print(names_df)
    names_df = names_df.reset_index(drop=True)
    print("names_df")
    print(names_df) #from here the nodes dataframe is empty
    names_df['node_id'] = names_df['node_id'].astype(str)
    name_mapping = dict(zip(names_df['node_id'], names_df['names']))

    # Füge die Namen zum index_genomes_df hinzu
    index_genomes_df['names'] = index_genomes_df['node_id'].map(name_mapping)
    print("index_genomes_df with names")
    print(index_genomes_df)



    index_genomes_df.to_csv(args.index_genomes, sep='\t', index=False, header=True) #'index_genomes_df.tsv'



    df_merged = nodes_df.merge(genomes_df, on='node_id', how='left') #merge nodes_df and genomes_df based on node_id
    print("df_merged")
    print(df_merged)
    #df_merged.to_csv('../evo_data/df_merged.tsv', sep='\t', index=False) #'df_merged.tsv'


    nodes_df['at_node'] = df_merged['at_node_y'].fillna(nodes_df['at_node']) #update column 4 in nodes_df with values from at_node from genomes_df
    #nodes_df.to_csv('../evo_data/nodes_df.tsv', sep='\t', index=False) #'nodes_df.tsv'


    def compute_sum_node(root, df):
        node_values = dict(zip(df["node_id"], df["at_node"]))
        sum_dict = {}
        range_dict = {}
        index = 0
        for node in postorder_iter(root):
            node_value = node_values.get(node.name, 0)
            child_sum = sum(sum_dict[child.name] for child in node.children)
            sum_dict[node.name] = node_value + child_sum
            if not node.children:
                range_dict[node.name] = (index, index + node_value)
            else:
                start = min(range_dict[child.name][0] for child in node.children)
                end = start + sum_dict[node.name]
                range_dict[node.name] = (start, end)
            index += node_value
        
        df["sum_node"] = df["node_id"].map(lambda x: sum_dict.get(x, 0))
        df["start_index"] = df["node_id"].map(lambda x: range_dict.get(x, (0, 0))[0])
        df["end_index"] = df["node_id"].map(lambda x: range_dict.get(x, (0, 0))[1])
        return df

    start = time.time()
    nodes_df = compute_sum_node(loadedtree,nodes_df) #call funktion on loaded tree and nodes_df

    end = time.time()
    print("time to traverse whole tree:", end - start)

    print("nodes_df")
    print(nodes_df)

    # table to tsv
    nodes_df.to_csv(args.summed_nodes, sep='\t', index=False) #'nodes_sum_table.tsv'


    #
    print(index_genomes_df)


if __name__ == "__main__":
    main()


