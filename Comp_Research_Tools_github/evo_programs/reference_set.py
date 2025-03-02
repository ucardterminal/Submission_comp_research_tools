# -*- coding: utf-8 -*-
import pandas as pd
import pickle
import random
import time
import numpy as np
from bigtree import Node, print_tree
from bigtree import Node, find_name
from bigtree import Node, tree_to_dot, postorder_iter
from bigtree import Node, findall, find_children ,add_dataframe_to_tree_by_path
from bigtree import Node, find, find_name, find_path, find_relative_path, find_full_path, find_attr
import argparse
import logging
import numpy as np

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='this program creates a reference set of genomes from the NCBI taxonomy tree')
    
    # Add arguments with default values
    parser.add_argument('--upperbound', type=str, default='131567',
                        help='provide a node_id for the upperbound, default it the cellular organisms=131567')
    parser.add_argument('--lowerbound', type=str, default="superkingdom",
                        help='provide a rank for the lowerbound, default is superkingdom')
    parser.add_argument('--num_of_seq', type=int, default=1,
                        help='state how many genomes you want to draw under each node of the lowerbound')
    parser.add_argument('--do_not_include', type=str, default=None,
                        help='state nodes, whose genomes should not be included in the reference set, you can provide multiple nodes separated by commas')
    parser.add_argument('--output_filename', type=str, default="../reference_set.tsv",
                        help='provide a filename for the reference set, default is reference_set.tsv')
    
    parser.add_argument(
        "--dryrun",
        action="store_true",
        help="Simulate the execution without making any changes."
    )
    
    # Parse the arguments
    args = parser.parse_args()
    dictionary = {}
    with open('../evo_data/results/ncbi_rank_dict.tsv', 'r') as file:
            lines = file.readlines()


            for line in lines:
 
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    number = int(parts[0])
                    string = parts[1]

                    dictionary[string] = number


    # Check if dryrun is enabled
    if args.dryrun:
        print("Only the size of the reference set with these parameters will be printed.")

        args_lowerbound =dictionary[args.lowerbound]
        args_lowerbound = int(args_lowerbound) #transform string to integer

        

        
        args = parser.parse_args()

        def parse_input(input_string):
            return [item.strip() for item in input_string.split(',')] #split the input string at the commas and remove the whitespaces
        if args.do_not_include is not None: #if the argument is provided
            do_not_include_list = parse_input(args.do_not_include)
        else:
            do_not_include_list = []


        arg_list = []

        for arg_name, arg_value in vars(args).items():
            # Check if the argument was explicitly provided or is using the default
            if arg_name in parser.parse_args([]).__dict__ and parser.parse_args([]).__dict__[arg_name] == arg_value:
                print(f"{arg_name}: {arg_value} (default)")
            else:
                print(f"{arg_name}: {arg_value}")
            
            # Add the argument to the list
            arg_list.append((arg_name, arg_value))



        with open("../evo_data/ncbitree.pickle", 'rb') as f:  #takes 3 seconds to load
            loadedtree = pickle.load(f)

        nodes_df = pd.read_csv('../evo_data/nodes_sum_table.tsv', delimiter='\t',header=0)
        nodes_df.columns = ['node_id', 'parent_id', 'rank','at_node','sum_node','start_index','end_index']

        nodes_df['node_id'] = nodes_df['node_id'].astype('int32')
        nodes_df['parent_id'] = nodes_df['parent_id'].astype('int32')
        nodes_df['rank'] = nodes_df['rank'].astype('int8')
        nodes_df['start_index'] = nodes_df['start_index'].astype('int32')
        nodes_df['end_index'] = nodes_df['end_index'].astype('int32')
        nodes_df['at_node'] = nodes_df['at_node'].astype('int32')
        nodes_df['sum_node'] = nodes_df['sum_node'].astype('int32')
    

        nodes_np = nodes_df.to_numpy()


        header = nodes_df.columns.to_numpy() #save the column names in a numpy array


        nodes_df['idx'] = nodes_df.loc[:, 'node_id'] # Set the node_id column  as index

        nodes_df.set_index('idx', inplace=True)



        index_genomes_sorted_df = pd.read_csv('../evo_data/index_genomes_df.tsv', delimiter='\t',header=0)
        index_genomes_sorted_df.columns = ['node_id','assembly_accession','names']




        post_order_list = []
        for children in postorder_iter(loadedtree):
            post_order_list.append(children.name) 
        post_order_list = list(map(int, post_order_list))

        index_genomes_df['node_id'] = pd.Categorical(index_genomes_df['node_id'], categories=post_order_list, ordered=True)
        index_refset_df = index_genomes_df.sort_values('node_id')
        neuer_index = pd.Index(range(len(index_refset_df)), name='index')
        index_refset_df = index_refset_df.reindex(neuer_index, copy=False)



        index_genomes_df['idx'] = index_genomes_df.loc[:, 'node_id']
        index_genomes_df.set_index('idx', inplace=True) #created a a new index for the dataframe with the node_ids as the indeces

        order_mapping = {node_id: i for i, node_id in enumerate(post_order_list)}


        # Then create a new column with the position values
        index_genomes_df['sort_key'] = index_genomes_df['node_id'].map(order_mapping)

        # Sort the DataFrame based on this new column
        index_genomes_df_sorted = index_genomes_df.sort_values('sort_key')
        index_genomes_df_sorted = index_genomes_df_sorted.drop('sort_key', axis=1)
        index_genomes_df_sorted = index_genomes_df_sorted.copy() 
        index_genomes_df_sorted.to_csv('../evo_data/index_genomes_df_sorted.tsv', sep='\t', index=False, header=True)




        Upper = find_name(loadedtree, str(args.upperbound))


        upperbound_descendants = []
        descendants = Upper.descendants
        for desc in descendants:
            upperbound_descendants.append(int(desc.name))



        filtered_upper_df = nodes_df[nodes_df['node_id'].isin(upperbound_descendants)]


        filtered_lower_upper_df = nodes_df[nodes_df['rank'] == args_lowerbound]

        filtered_lower_upper_df.to_csv('filtered_lower_upper_df.tsv', sep='\t', index=False, header=True)


                
        def visit_children(node):
            list_check = []
            for children in node.children:
                #print(children.node_name)
                try:
                    node_value = int(children.node_name) # gets the children and typcasts them to int

                    rank_idx = nodes_df.loc[node_value, 'rank']

                    if children.node_name in do_not_include_list:

                        continue
                    if rank_idx >= args_lowerbound or len(children.children) == 0:
                        list_check.append(children.node_name)

                        continue  # skips the resursive call for this node
                except ValueError:
                    print("error")
                    # if node_name is not a number, then
                    pass
                    
                
                child_results = visit_children(children) #only recursively call if the node does not meet the condition
                list_check.extend(child_results) #extend the list with the results from the recursive call
                
            return list_check 

        start = time.time()
        startnode = find_name(loadedtree, str(args.upperbound))  # has to be type casted to string eventhough it is initialized as a sting
        lowerbound_set = visit_children(startnode)
        end = time.time()


        leaves = loadedtree.leaves
        leaveslist = []
        for leaf in leaves:
            leaveslist.append(leaf.name)

        #print(do_not_include_list)
        main_list = list(set(leaveslist) - set(lowerbound_set)) #get the difference between the two lists


        lowerbound_set_int = list(map(int, lowerbound_set)) #typecast all the elements in the list to int

        lowerboundset_df = nodes_df.loc[nodes_df.index.isin(lowerbound_set_int)]

        lowerboundset_df =lowerboundset_df.copy() #to avaid the warning it has to be an actual dataframe not a slice
        lowerboundset_df.to_csv('../evo_data/lowerbound_set_df.tsv', sep='\t', index=False, header=True)

        def generate_random_numbers(x, y, n):
            if y - x <= 0:
                return []  # Return an empty list if y - x = 0
            elif y - x == 1:
                return [x] #if difference between x and y is 1, return x
            count = min(n, y - x)   
            return random.sample(range(x, y + 1), count) #with random.sample, duplicates are not allowed



        lowerboundset_df.loc[:,'final_index'] = lowerboundset_df.apply(lambda row: generate_random_numbers(row['start_index'], row['end_index'], args.num_of_seq), axis=1)

        combined_list = sum(lowerboundset_df['final_index'].tolist(), [])

        print("num of genomes for the referenceset with these settings : ",len(combined_list))
        # Simulate the program's actions here
    else: #if not dryrun

        args_lowerbound =dictionary[args.lowerbound]
        print(type(args_lowerbound))
        #args_lowerbound = int(args_lowerbound) #transform string to integer


        # Parse the arguments
        args = parser.parse_args()

        def parse_input(input_string):
            return [item.strip() for item in input_string.split(',')]
        if args.do_not_include is not None: #if the argument is provided
            do_not_include_list = parse_input(args.do_not_include)
        else:
            do_not_include_list = []


        arg_list = []

        for arg_name, arg_value in vars(args).items():
            if arg_name in parser.parse_args([]).__dict__ and parser.parse_args([]).__dict__[arg_name] == arg_value:
                print(f"{arg_name}: {arg_value} (default)")
            else:
                print(f"{arg_name}: {arg_value}")
            
            # Add the argument to the list
            arg_list.append((arg_name, arg_value))



        with open("../evo_data/ncbitree.pickle", 'rb') as f:  #takes 3 seconds to load
            loadedtree = pickle.load(f)

        nodes_df = pd.read_csv('../evo_data/nodes_sum_table.tsv', delimiter='\t',header=0)
        nodes_df.columns = ['node_id', 'parent_id', 'rank','at_node','sum_node','start_index','end_index']

        nodes_df['node_id'] = nodes_df['node_id'].astype('int32')
        nodes_df['parent_id'] = nodes_df['parent_id'].astype('int32')
        nodes_df['rank'] = nodes_df['rank'].astype('int8')
        nodes_df['start_index'] = nodes_df['start_index'].astype('int32')
        nodes_df['end_index'] = nodes_df['end_index'].astype('int32')
        nodes_df['at_node'] = nodes_df['at_node'].astype('int32')
        nodes_df['sum_node'] = nodes_df['sum_node'].astype('int32')

        nodes_df.to_csv('nodes_with_sum_node_df.tsv', sep='\t', index=False, header=True)
    

        nodes_np = nodes_df.to_numpy()
        print("nodes_np")
        print(nodes_np)

        # Header (Spaltennamen) speichern
        header = nodes_df.columns.to_numpy()
        print("Header:")
        print(header)


        # Setze die node_id-Spalte als Index
        nodes_df['idx'] = nodes_df.loc[:, 'node_id']

        nodes_df.set_index('idx', inplace=True)



        print(nodes_df)
        index_genomes_df = pd.read_csv('../evo_data/index_genomes_df.tsv', delimiter='\t',header=0)
        index_genomes_df.columns = ['node_id','assembly_accession','names']

        print("index_genomes_df")
        print(index_genomes_df)



        post_order_list = []
        for children in postorder_iter(loadedtree):
            post_order_list.append(children.name) 
        post_order_list = list(map(int, post_order_list))

        index_genomes_df['node_id'] = pd.Categorical(index_genomes_df['node_id'], categories=post_order_list, ordered=True)
        index_refset_df = index_genomes_df.sort_values('node_id')
        neuer_index = pd.Index(range(len(index_refset_df)), name='index')
        index_refset_df = index_refset_df.reindex(neuer_index, copy=False)

        print("index_refset_df")
        print(index_refset_df)

        index_genomes_df['idx'] = index_genomes_df.loc[:, 'node_id']
        index_genomes_df.set_index('idx', inplace=True) 

        order_mapping = {node_id: i for i, node_id in enumerate(post_order_list)}



        # Then create a new column with the position values
        index_genomes_df['sort_key'] = index_genomes_df['node_id'].map(order_mapping)
        print("index_genomes_df")
        print(index_genomes_df)
        # Sort the DataFrame based on this new column
        index_genomes_df_sorted = index_genomes_df.sort_values('sort_key')
        # Optionally, drop the temporary sort_key column if you don't need it anymore
        index_genomes_df_sorted = index_genomes_df_sorted.drop('sort_key', axis=1)
        index_genomes_df_sorted = index_genomes_df_sorted.copy()
        index_genomes_df_sorted.to_csv('../evo_programs/index_genomes_df_sortedupdate.tsv', sep='\t', index=False, header=True)



        Upper = find_name(loadedtree, str(args.upperbound))
        print("Upper")
        print(Upper)

        upperbound_descendants = []
        descendants = Upper.descendants
        for desc in descendants:
            upperbound_descendants.append(int(desc.name))



        filtered_upper_df = nodes_df[nodes_df['node_id'].isin(upperbound_descendants)]
        print(filtered_upper_df)

        filtered_lower_upper_df = nodes_df[nodes_df['rank'] == args_lowerbound]
        print("filtered_lower_upper_df")
        print(filtered_lower_upper_df)
        print(len(filtered_lower_upper_df))
        filtered_lower_upper_df.to_csv('filtered_lower_upper_df.tsv', sep='\t', index=False, header=True)


                
        def visit_children(node):
            list_check = []

            for children in node.children:
                #print(children.node_name)
                try:
                    node_value = int(children.node_name) # gets the children and typcasts them to int
                    rank_idx = nodes_df.loc[node_value, 'rank']
                    if children.node_name in do_not_include_list:

                        continue
                    if rank_idx >= args_lowerbound or len(children.children) == 0:
                        list_check.append(children.node_name)
                        
                        continue  # skips the resursive call for this node
                except ValueError:
                    print("error")
                    # if node_name is not a number, then
                    pass
                    
                
                child_results = visit_children(children) #only recursively call if the node does not meet the condition
                list_check.extend(child_results) #extend the list with the results from the recursive call
                
            return list_check 

        start = time.time()
        startnode = find_name(loadedtree, str(args.upperbound))  # has to be type casted to string eventhough it is initialized as a sting
        lowerbound_set = visit_children(startnode)
        end = time.time()

        #print(result)
        print(len(lowerbound_set))
        print("time taken :" ,end-start)

        leaves = loadedtree.leaves
        leaveslist = []
        for leaf in leaves:
            leaveslist.append(leaf.name)
        print(len(leaveslist))

        #print(do_not_include_list)
        main_list = list(set(leaveslist) - set(lowerbound_set)) #get the difference between the two lists
        print("len of main list")
        print(len(main_list))
        print("main list")
        #print(main_list)

        lowerbound_set_int = list(map(int, lowerbound_set)) #typecast all the elements in the list to int

        lowerboundset_df = nodes_df.loc[nodes_df.index.isin(lowerbound_set_int)]
        print(lowerboundset_df)
        lowerboundset_df =lowerboundset_df.copy() #to avaid the warning it has to be an actual dataframe not a slice

        def generate_random_numbers(x, y, n):
            if y - x <= 0:
                return []  # Return an empty list if y - x = 0
            elif y - x == 1:
                return [x] #if difference between x and y is 1, return x
            count = min(n, y - x)   
            return random.sample(range(x, y + 1), count) #with random.sample, duplicates are not allowed


        print(lowerboundset_df.dtypes)

        lowerboundset_df.loc[:,'final_index'] = lowerboundset_df.apply(lambda row: generate_random_numbers(row['start_index'], row['end_index'], args.num_of_seq), axis=1)
        print("lowerboundset_df")
        print(lowerboundset_df)

        combined_list = sum(lowerboundset_df['final_index'].tolist(), []) #extracts al the list values of the final_index column to a list

        print("num of genomes for the referenceset with these settings : ",len(combined_list))
    
        #typecast the combined list to int
        combined_list = list(map(int, combined_list))
        print("combined_list")
        # combined_list pandas dataframe
        randomly_chosen_represenatives_df = pd.DataFrame(combined_list, columns=['node_id'])
        randomly_chosen_represenatives_df.to_csv('../evo_data/randomly_chosen_represenatives_df.tsv', sep='\t', index=False, header=True)

        #if the element of the list are in the column node_id of the dataframe index_genomes_df_sorted, then return the row
        #reindex index_genomes_df_sorted based on the linear index
        index_genomes_df_sorted = index_genomes_df_sorted.reset_index(drop=True)
        refset_assemblys_df = index_genomes_df_sorted .loc[combined_list] #####
        refset_assemblys_df.to_csv('refset_assemblystest_df.tsv', sep='\t', index=False, header=True)


        #refset_assemblys_df = index_genomes_df_sorted.loc[index_genomes_df_sorted.index.isin(combined_list)]
        print("dataframe of the final reference set")
        print(refset_assemblys_df)
        refset_assemblys_df = refset_assemblys_df[['names', 'node_id', 'assembly_accession']] #change the appearance of the columns in the dataframe
        refset_assemblys_df.to_csv(args.output_filename, sep='\t', index=False, header=True) #save the dataframe to tsv


if __name__ == "__main__":
    main()