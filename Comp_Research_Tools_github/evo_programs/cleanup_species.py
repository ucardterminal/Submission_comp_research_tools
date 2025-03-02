# -*- coding: utf-8 -*-
import pandas as pd
import argparse
def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='clean up the UTF- error if there is any')
    
    # Add arguments with default values
    parser.add_argument("--file_to_be_cleaned", default='../evo_data/species.txt',help="filename of input file")
    parser.add_argument("--cleaned_file", default='../evo_data/results/species_cleaned.txt',help="filename of output file")

    args = parser.parse_args()


    def check_utf8(file_path):
        error_lines = []
        
        with open(file_path, "rb") as f:  # binary mode to check individual bytes
            for line_number, line in enumerate(f, start=1):
                try:
                    line.decode("utf-8")  # try to decode as UTF-8
                except UnicodeDecodeError as e:
                    error_lines.append((line_number, e.reason))  # save line number + error description in the list
        
        if error_lines:
            print(f"There were {len(error_lines)} faulty lines:")
            for line_number, reason in error_lines:
                print(f"line {line_number}: {reason}")
        else:
            print("All characters are in UTF-8 format") 


    check_utf8(args.file_to_be_cleaned)

    with open(args.file_to_be_cleaned, "r", encoding="utf-8", errors="replace") as f: #file_to_be_cleaned   "species.txt"
        cleaned_content = f.read()

    with open(args.cleaned_file, "w", encoding="utf-8") as cleaned_file: #cleaned_file   "species_cleaned.txt"
        cleaned_file.write(cleaned_content)

    print("cleaned file saved as species_cleaned.txt 'species_cleaned.txt' ")


if __name__ == "__main__":
    main()
