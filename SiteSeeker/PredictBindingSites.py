#!/usr/bin/env python

from joblib import load
import argparse
import os
import sys
from Bio.PDB.PDBIO import PDBIO
import SiteSeeker.residue_functions as residue_functions
from sklearn.preprocessing import OneHotEncoder
import pandas as pd
import numpy as np
from sklearn import svm
from sklearn.metrics import accuracy_score, precision_score, recall_score
import pkg_resources

# Loading the model and the scaler:
# Get the directory of this module
module_dir = os.path.dirname(__file__)

# Construct the path to the model.joblib file relative to the module directory
model_path = os.path.join(module_dir, 'model.joblib')
# Load the model
model = load(model_path)

# Construct the path to the model.joblib file relative to the module directory
scaler_path = os.path.join(module_dir, 'scaler.pkl')
# Load the scaler
scaling_params = load(scaler_path)


def ensure_columns(df):
    """Ensures that all columns from this specific list of column names are present in the inputted DataFrame 
    and makes sure the columns are in the correct order for the machine learning model.
    If any columns are missing, they are created with all zero values."""
    # List of column names that should be present in the DataFrame in the correct order
    column_names = ['kd_hydrophobicity', 'rose_hydrophobicity', 'charge', 'polarizability',
                    'molar_mass', 'volume_vdw', 'pI', 'steric_parameter', 'ASA', 'phi',
                    'psi', 'NH_O_1_energy', 'O_NH_1_energy', 'NH_O_2_energy',
                    'O_NH_2_energy', 'neighbor_density', 'ASA_Zscore', 'residue_A',
                    'residue_C', 'residue_D', 'residue_E', 'residue_F', 'residue_G',
                    'residue_H', 'residue_I', 'residue_K', 'residue_L', 'residue_M',
                    'residue_N', 'residue_P', 'residue_Q', 'residue_R', 'residue_S',
                    'residue_T', 'residue_V', 'residue_W', 'residue_Y',
                    'secondary_structure_-', 'secondary_structure_B',
                    'secondary_structure_E', 'secondary_structure_G',
                    'secondary_structure_H', 'secondary_structure_I',
                    'secondary_structure_P', 'secondary_structure_S',
                    'secondary_structure_T']
    # Get any missing columns
    missing_columns = set(column_names) - set(df.columns)
    # If any columns are missing, create them with all zero values
    for col in missing_columns:
        df[col] = 0
    # Reorder the columns to match the expected order
    df = df.reindex(columns=column_names)
    return df

def groupBindingSites(featuresTable, predictions, max_dist = 4):
    """Groups the binding residues into binding sites based on the distance between them."""
    # Get the binding residues (binding_prediction == 1) and concatenate the predictions to the features table
    BindingResidues = pd.concat([featuresTable, predictions], axis=1)[predictions["binding_prediction"] == 1]
    # List of lists to store the binding sites
    BindingSites = []
    # Iterate over the binding residues
    for i in range(0, len(BindingResidues)): 
        # Get the residue instance
        residue1 = BindingResidues.iloc[i]['residue_instance']
        if i == 0:
            # If it is the first residue, create a new binding site with the residue in it
            BindingSites.append([residue1])
        else:
            # If it's not the first residue, check if it is close to any of the existing binding sites
            sitesList = []
            for site_j in range(0, len(BindingSites)):
                binds = False
                # If the residue is close to any of the residues in the binding site, it's considered to be close to the site
                for residue2 in BindingSites[site_j]:
                    dist = residue_functions.minimum_distance(residue1, residue2)
                    if dist <= max_dist:
                        binds = True
                        break
                if binds:
                    # Add all sites it is close with to this list
                    sitesList.append(site_j)
            # If the residue is close to only one site, add it to that site
            if len(sitesList) == 1:
                BindingSites[sitesList[0]].append(residue1)
            # If the residue is not close to any site, create a new site with it
            elif len(sitesList) == 0:
                BindingSites.append([residue1])
            # If the residue is close to more than one site, merge the sites
            else:
                newSite = [residue1]
                for index in sitesList[::-1]:
                    newSite.extend(BindingSites[index])
                    # Remove the old unmerged sites from the list
                    BindingSites.pop(index)
                # Add the new merged site
                BindingSites.append(newSite)
    return(BindingSites)


def SiteSeeker(input_file, known_binding = False, verbose = False, output = None, max_distance = 4):
    """
    Predicts binding sites in a protein structure (from a PDB file) based on a SVM model which considers a set of 
    structural and physicochemical features, trained on a dataset of protein structures with known binding sites.
    Outputs a PDB file with the binding residues and a text file with the binding sites and the residues involved
    in them.
    If known_binding is set to True, caluculates the accuracy, precision and recall based on the residues that bind
    to the heteroatoms in the input file (excluding water molecules).
    """
    # Getting the PDB code from the input file name (if the file name is not code.pdb, the code will be whatever comes before .pdb)
    pdb_code = input_file.split('/')[-1].split('.pdb')[0]
    # Extracting the features from the input file to get the features table
    if verbose:
        print("Extracting features from input file...")
    features = residue_functions.residue_features(pdb_code, input_file, binding=known_binding)
    if verbose:
        print("Predicting binding sites...")
    # Encoding categorical features into binary
    one_hot_encoder = OneHotEncoder()
    encoded_features = one_hot_encoder.fit_transform(features[['residue', 'secondary_structure']])
    # Converting the encoded features into a DataFrame
    encoded_df = pd.DataFrame(encoded_features.toarray(), columns=one_hot_encoder.get_feature_names_out(['residue', 'secondary_structure']))
    # Concatenating the encoded features with the original DataFrame
    dataset_encoded = pd.concat([features, encoded_df], axis=1)
    # Removing the original categorical features
    dataset_encoded.drop('residue', axis=1, inplace=True)
    dataset_encoded.drop('secondary_structure', axis=1, inplace=True)
    # Getting the features for the machine learning model
    if known_binding:
        # If known_binding is True, extract also y (binding / non-binding)
        x_features = dataset_encoded.iloc[:, 5:]
        y_features = dataset_encoded.iloc[:, 4]
    else:
        # If not, just get the features to predict
        x_features = dataset_encoded.iloc[:, 4:]
    # Making sure the table has all the columns needed for the model and are in the right order
    x_features_complete = ensure_columns(x_features) 
    # Scaling the features using the saved scaler
    x_features_scaled = scaling_params.transform(x_features_complete)
    # Converting the scaled features back to a DataFrame
    x_features_scaled_df = pd.DataFrame(x_features_scaled, columns=x_features_complete.columns, index=x_features_complete.index)
    # Predicting the binding sites
    prediction = model.predict(x_features_scaled_df)
    if verbose:
        print("Obtaining results...")
    # Grouping the binding residues into binding sites
    BindingSites = groupBindingSites(features, pd.DataFrame(prediction, columns=['binding_prediction']), max_distance)
    # Determine output file names
    if output is None:
        # If no output directory is specified, save the files in the current directory as pdb_code_SiteSeeker.pdb and pdb_code_SiteSeeker.txt
        outputpdb = "./" + pdb_code + '_SiteSeeker.pdb'
        outputtxt = "./" + pdb_code + '_SiteSeeker.txt'
    else:
        # Else, save the files in that directory with the specified name
        outputpdb = output + '.pdb'
        outputtxt = output + '.txt'
    # Save the binding residues to a PDB file
    with open(outputpdb, 'w') as output_file:
        # It will be a PDB file containing only the residues in binding sites (to be opened alongside the original PDB file in Chimera)
        for num in range(0, len(prediction)):
            if prediction[num] == 1:
                io = PDBIO()
                io.set_structure(features.loc[num, 'residue_instance'])
                io.save(output_file)
    # Save the binding sites to a text file
    with open(outputtxt, 'w') as output_file:
        for num in range(0, len(BindingSites)):
            # Save each binding site 
            print("Binding Site " + str(num+1), file= output_file)
            # Sort the residues by chain and residue number
            BindingSites[num].sort(key = lambda x: (x.get_parent().get_id(), x.get_id()[1]))
            # Save each residue: chain, residue, number
            for residue in BindingSites[num]:
                residueID = residue.get_id()
                residuename = residue.get_resname()
                chain = residue.get_parent().get_id()
                print(chain, residuename, residueID[1], sep = "\t", file = output_file)
    print(str(len(BindingSites)) + ' binding sites were predicted. Results can be found in ' + outputtxt + ', and the binding residues can be visualized with the file ' + outputpdb + ' using Chimera.')
    # If known_binding is True, calculate the accuracy, precision and recall and print them
    if known_binding:
        accuracy = accuracy_score(y_features, prediction)
        precision = precision_score(y_features, prediction)
        recall = recall_score(y_features, prediction)
        print("Accuracy: " + str(round(accuracy,4)))
        print("Precision: " + str(round(precision,4)))
        print("Recall: " + str(round(recall,4)))
    # Return the prediction in case the function is called from another script
    return(pd.DataFrame(prediction, columns=['binding_prediction']))

# Main function to be executed when the script is run (saved in a function to be installed with setuptools)
def main():
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="""
                                     Predicts binding sites in a protein structure (from a PDB file) based on a SVM model which considers a set of 
                                     structural and physicochemical features, trained on a dataset of protein structures with known binding sites.
                                     Outputs a PDB file with the binding residues and a text file with the binding sites and the residues involved
                                     in them.
                                     If knownbind is set, caluculates the accuracy, precision and recall based on the residues that bind
                                     to the heteroatoms in the input file (excluding water molecules).
                                     """)
    # Add input arguments
    parser.add_argument('--input', '-i', type=str, help='Input file in PBD format (with extension .pdb)', required=True)
    parser.add_argument('--output', '-o', type=str, help='Desired path and base filename of the output files (without the file extension)', required=False)
    parser.add_argument('--verbose', '-v', action='store_true', help='Prints additional information during the execution', required=False)
    parser.add_argument('--maxdist', '-d', type=int, help='Maximum distance between the closest residues in two binding sites for these to be considered separate binding sites (default: 4 Å)', required=False)
    parser.add_argument('--knownbind', '-k', action="store_true", help='The input PDB file contains the interaction with the ligand (with the flag HETATM) and I would like to know the accuracy of SiteSeeker\'s prediction',  required=False)
    # Parse the command-line arguments
    args = parser.parse_args()
    # Access the input argument
    input_arg = args.input
    if input_arg:
        if os.path.exists(input_arg):
            if os.path.isfile(input_arg) and input_arg.endswith(".pdb"):
                if args.verbose:
                    print("Reading info from ", input_arg)
            else:
                print(f"The file {input_arg} must be a PDB file in .pdb format.", file=sys.stderr)
                sys.exit(1)
        else:
            print(f"The file {input_arg} does not exist.", file=sys.stderr)
            sys.exit(1)
    else:
        print("The required argument '--input' is missing.", file=sys.stderr)
        sys.exit(1)
    # Save the maximum distance between residues in a different binding site if given; default to 4 A otherwise
    if args.maxdist:
        max_distance = args.maxdist
    else:
        max_distance = 4
    # Save the output path if given; default to None otherwise
    if args.output: 
        output_dir = args.output
    else:
        output_dir = None
    # Call the SiteSeeker function with the input arguments
    SiteSeeker(input_arg, known_binding= args.knownbind, verbose = args.verbose, output= output_dir, max_distance = max_distance)

# If the script is run directly, call the main function
if __name__ == '__main__':
    main()

    

