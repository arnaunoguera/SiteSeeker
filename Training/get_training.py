from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBList
from Bio.PDB.DSSP import DSSP
import numpy as np
import pandas as pd
import aminoacid_info as aa
import residue_functions as functions
import os

cont  = 0
first = True
output_file = 'training_file.tsv'

# Iterate through the files of a folder
folder_path = './pdb_files/zip_Arnau/'
for filename in os.listdir(folder_path):
    file_path = os.path.join(folder_path, filename)
    #Check if it is a file
    if os.path.isfile(file_path) and file_path.endswith('.pdb'):
        # Process the file here
        pdb_code = file_path.split('/')[-1].replace('.pdb', '')
        print("Trying " + pdb_code)
        try:
            if first:
                results_table = functions.residue_features(pdb_code, file_path)
                results_table.to_csv(output_file, sep='\t', index=False)
                first = False
            else: 
                new_data = functions.residue_features(pdb_code, file_path)
                with open(output_file, 'a') as f:
                    new_data.to_csv(f, sep='\t', header=False, index=False)
            print(str(cont) + ".File "  + file_path.split('/')[-1] + " succesfully read")
        except Exception as e:
            print(f"Error processing {filename}: {e}")
        finally:
            cont += 1
       
print("All files processed and saved to ", output_file)