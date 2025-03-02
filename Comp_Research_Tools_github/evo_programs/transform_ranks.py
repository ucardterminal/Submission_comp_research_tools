# -*- coding: utf-8 -*-
import argparse


def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='this program transforms the nodes.dmp file, it replaces the ranks_stings with their keys from the ranks_from_ncbi.tsv file')
    
    # Add arguments with default values
    parser = argparse.ArgumentParser()
    parser.add_argument("--input1", default='../evo_data/ranks_from_ncbi.tsv',help="filename of input file")
    parser.add_argument("--input2", default='../evo_data/nodes.dmp',help="filename of input file")
    parser.add_argument("--output1",default='../evo_data/results/ncbi_rank_dict.tsv',  help="filename of nodes output file")
    parser.add_argument("--output2",default='../evo_data/results/nodes_copy.tsv',  help="filename of dict output file")
    args = parser.parse_args()




    def create_tree_dictionary(filename):
        ranks = {}
        with open(filename, 'r') as file:
            for index, line in enumerate(file): #for each line
                string_part = line.strip().split()[0]# plits the line and take only the first element
                cleaned_string = string_part.replace('_', ' ') #Replace underscores with spaces
                # Add to dictionary with row number as key
                cleaned_string = cleaned_string.lower()
                ranks[index] = cleaned_string
        return ranks


    # Using the function the function
 
    ranks = create_tree_dictionary(args.input1) #default='ranks_from_ncbi.tsv'

    # Print the resulting dictionary (optional)
    for key, value in ranks.items():
        print(f"{key}: {value}")


    with open(args.output1, "w") as file: #default='ncbi_rank_dict.tsv'
        for key, value in ranks.items():
            file.write(f"{key}\t{value}\n")



    with open(args.input2, "r") as infile, open(args.output2, "w") as outfile: #default input2 : nodes.dmp output2: "nodes_copy.tsv"
        for line in infile:
            fields = [field.strip() for field in line.split('|')] # Split the line by '|'
            if len(fields) >= 3:# Check if there are at least 3 fields
                rank = fields[2]#get the rank
                if rank in ranks.values():# Check if the rank is in the dictionary
                    #replace the rank with the corresponding key
                    rank_key = next(key for key, value in ranks.items() if value == rank)
                    fields[2] = str(rank_key)
                    
                outfile.write("\t".join(fields[:3]) + "\n") # Write the first 3 fields separated by tabs



if __name__ == "__main__":
    main()