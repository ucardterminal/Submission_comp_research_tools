# -*- coding: utf-8 -*-
import pandas as pd
import time
from bigtree import dataframe_to_tree_by_relation
import pickle
start = time.time()
import argparse

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='this program creates a tree from a dataframe')
    
    # Add arguments with default values
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default='../evo_data/nodes_copy.tsv',help="filename of input file")
    parser.add_argument("--output", default='../evo_data/results/ncbitree.pickle',help="filename of output file")
    args = parser.parse_args()


    df = pd.read_csv(args.input, delimiter='\t',header=None)
    df.columns = ['node_id', 'parent_id', 'rank']


    child_df = pd.DataFrame(0, index=range(0), columns=range(3))
    child_df.columns = ['node_id', 'parent_id', 'rank']
    child_df.iloc[:] = -1
    filtered_rows = df[(df.iloc[:, 1] == 1)| (df.iloc[:, 1] == 2787823)| (df.iloc[:, 1] == 2787854)| (df.iloc[:, 1] == 10239)| (df.iloc[:, 1] == 131567)]



    # Die gefilterten Zeilen zu df2 hinzufügen
    child_df = pd.concat([child_df, filtered_rows], ignore_index=True)
    #hier wird filtered rows an new_df angehängt

    print(child_df)


    sorted_child_df = child_df.sort_values(by=[child_df.columns[1],child_df.columns[0]])
    print(sorted_child_df.head(87))
    sorted_child_df['path'] = -1

    print("sorted_child_df")
    print(sorted_child_df)


    def find_paths_and_update(df):
        for index, row in df.iterrows():
            current_node = row['node_id']
            path = []  # Temporäre Liste für den aktuellen Pfad

            while current_node != 1:  # Solange der aktuelle Knoten nicht 1 ist
                path.append(int(current_node))  # Füge den aktuellen Knoten zur Liste hinzu und konvertiere zu int
                # Suche den parent_id des aktuellen Knotens
                parent_row = df[df['node_id'] == current_node]
                if not parent_row.empty:
                    current_node = parent_row['parent_id'].values[0]  # Aktualisiere den aktuellen Knoten
                else:
                    break  # Breche ab, wenn kein Elternknoten gefunden wird

            # Füge die Wurzel (1) hinzu und erstelle den Pfad im gewünschten Format

            end=time.time()
            print(end-start)

            path.append(1)  # Füge die Wurzel (1) hinzu
            path_str = '/'.join(map(str, reversed(path)))  # Erstelle den Pfad-String im Format 1/278785/28384/
            
            # Aktualisiere die 'path'-Spalte im DataFrame
            df.at[index, 'path'] = path_str  # Setze den Pfad-String in der entsprechenden Zeile

        return df

    find_paths_and_update(sorted_child_df)
    print("sorted_child_df")
    print(sorted_child_df)
    sorted_child_df.iloc[:, [0, 1]] = sorted_child_df.iloc[:, [0, 1]].astype(str)
    sorted_child_df.iloc[0, 1] = None
    print("sorted_child_df")
    print(sorted_child_df)

    root = dataframe_to_tree_by_relation(sorted_child_df)
    #root.show(attr_list=["node_id"])


    print(df)
    neuer_df = df.iloc[:10000].copy()
    print("neuer_df")
    print(neuer_df) 
    cols = df.columns[:2]
    df[cols] = df[cols].astype(str)
    #neuer_df.iloc[0, 1] = None
    print("aktualisiert")
    print(neuer_df)


    largerchild_df = df[df.iloc[:, 1].isin(sorted_child_df.iloc[:, 0])] #fügt df-neu alle zeilen von df hinzu, die in column 1 sortedchild_df sind
    print("largerchild_df")
    print(largerchild_df)



    #root = dataframe_to_tree_by_relation(sorted_child_df)
    largerchild_df[cols] = largerchild_df[cols].astype(str)
    largerchild_df.iloc[0, 1] = None



    df[cols] = df[cols].astype(str)
    df.iloc[0, 1] = None
    root = dataframe_to_tree_by_relation(largerchild_df)
    #root.show(attr_list=["node_id"])


    with open(args.output, "wb") as f:
        pickle.dump(root, f)

    print("largerchild_df")
    print(largerchild_df)
    print("done")



if __name__ == "__main__":
    main()