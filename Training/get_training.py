#Importing the necessary modules.

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBList
from Bio.PDB.DSSP import DSSP
import numpy as np
import pandas as pd
import aminoacid_info as aa
import residue_functions as functions
import os

# Defining the variables.

cont  = 1
first = True
output_file = 'training_file.tsv'

# Iterating through the files of a folder and checking if they are a pdb file.
folder_path = './pdb_files/zip_Arnau/'
for filename in os.listdir(folder_path):
    file_path = os.path.join(folder_path, filename)
    if os.path.isfile(file_path) and file_path.endswith('.pdb'):

        # Processing the filename to obtain the pdb code of the file.

        pdb_code = file_path.split('/')[-1].replace('.pdb', '')
        print("Trying " + pdb_code)

        #Since there are some pdb files that don't work in dssp, using try and except to skip the files that give errors.

        try:

            #If it is the first file, calculating the residue features of the protein and creating a csv file with said features.

            if first:
                results_table = functions.residue_features(pdb_code, file_path)
                results_table.to_csv(output_file, sep='\t', index=False)
                first = False

            #If it is not the first file, calculating the residue features of the protein and apending them to the previous csv file.
                
            else: 
                new_data = functions.residue_features(pdb_code, file_path)
                with open(output_file, 'a') as f:
                    new_data.to_csv(f, sep='\t', header=False, index=False)
            print(str(cont) + ".File "  + file_path.split('/')[-1] + " succesfully read")

        #Skiping the files that give errors during text processing (tipically problems with dssp)

        except Exception as e:
            print(f"Error processing {filename}: {e}")

        #Adding 1 to a cont to keep a track of the number of already processed files independently of if they gave an error or not. 

        finally:
            cont += 1
print("All files processed and saved to ", output_file)