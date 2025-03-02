# -*- coding: utf-8 -*-
import pandas as pd
import time
from bigtree import dataframe_to_tree_by_relation
import pickle
import argparse


def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='this program creates a tree from a dataframe')
    
    # Add arguments with default values
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default='../evo_data/nodes_copy.tsv',help="filename of input file")
    parser.add_argument("--output", default='../evo_data/results/ncbitree.pickle',help="filename of output file")
    args = parser.parse_args()

    start = time.time()

    df = pd.read_csv(args.input, delimiter='\t',header=None) #default='nodes_copy.tsv'
    df.columns = ['node_id', 'parent_id', 'rank']

    cols = df.columns[:2]
    df[cols] = df[cols]#.astype(str)

    copy_df = df.sort_values(by=[df.columns[1],df.columns[0]])
    df = copy_df.copy()
    print("dataframe build")

    df[cols] = df[cols].astype(str)
    df.iloc[0, 1] = None

    end = time.time()
    print("tree is being built")
    print(end-start)

    root = dataframe_to_tree_by_relation(df) #the tree is being built here
    #root.show(attr_list=["node_id"]) #hier wird der baum angezeigt 

    end = time.time()
    print("tree build")
    print(end-start)

    with open("ncbitree.pickle", "wb") as f: #the tree is being saved here
        pickle.dump(root, f)

    end = time.time()
    print("saves as pickle",end-start)

    print("df")
    print(df)
    end = time.time()
    print("final time")
    print(end-start)
    print("done")

if __name__ == "__main__":
    main()