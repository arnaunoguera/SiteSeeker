#Importing the necessary modules.

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBList
from Bio.PDB.DSSP import DSSP
import numpy as np
import pandas as pd
import SiteSeeker.aminoacid_info as aa

#Function to get the distance between two atoms.

def get_distance(atom1, atom2):
    
    #Getting the coordinates of the two atoms and calculating the distance between them.

    coord1 = np.array(atom1.get_coord())
    coord2 = np.array(atom2.get_coord())
    dist = np.sqrt(np.sum((coord1-coord2)**2))

    #Returning the distance.

    return dist

#Function to get the minimum distance between two residues.

def minimum_distance(residue1, residue2):
    first = True

    #Iterating trough the atoms of the two residues.
    for res1atom in residue1.get_iterator():
        for res2atom in residue2.get_iterator():

            #Calculating the distance between a pair atoms.

            dist = get_distance(res1atom,res2atom)

            #If it is the first pair of atoms or the distance is lower than the minimum distance, storing the distance in the minimum distance.
            if first == True or dist < min_dist:
                first = False
                min_dist = dist
    
    #Returning the minimum distance.
                
    return min_dist

#Function to calculate if a residue is concave or convex with respect to the ones surrounding it, as well as neighbor density. 

def calculate_density_and_convexity(structure, reference_residue, search_radius, ASA, dssp):

    #Defining the necessary variables. 

    neighbor_count = 0
    neighbor_list_ASA = []

    #Iterating trough the residues of a given structure and skipping heteroatoms.

    for residue2 in structure.get_residues():
        residue2ID = residue2.get_id() 
        if residue2ID[0][0] == ' ':
            
            #Calculating the minimum distance between the residue and the reference residue.

            distance = minimum_distance(reference_residue, residue2)

            #If the distance is equal or lower than a given radious, the residue is considered a neighbor of the reference residue.

            if distance <= search_radius:
                neighbor_count += 1

                #Obtaining  the chainID and searching the residue 2 in the dssp table. Storing the ASA in a list of ASAs.

                chainID = residue2.get_parent().get_id()
                dssp_key = (chainID, residue2ID)
                try:
                    residue2ASA = dssp[dssp_key][3]
                except KeyError:
                    neighbor_count -= 1
                    continue
                neighbor_list_ASA.append(residue2ASA)

    # Calculating volume of a sphere with a radius equal to the search radious.
                
    sphere_volume = (4/3) * np.pi * (search_radius ** 3)
    
    # Calculating neighbor density.

    neighbor_density = neighbor_count / sphere_volume

    #Calculating ASA_Z.

    ASA_Z = ((ASA-np.mean(neighbor_list_ASA))/np.std(neighbor_list_ASA))

    #Returning neighbour density and the ASA_Z.
    
    return(neighbor_density, ASA_Z)

#Function to calculate the features of a given residue.

def residue_features(pdb_code, pdb_file, binding=True):

    #Starting PDBParser

    p = PDBParser(PERMISSIVE=1, QUIET=1)
    
    #Obtaining the protein structure of a given file.

    structure = p.get_structure(pdb_code, pdb_file)

    #Defining the model and using dssp to calculate secondary structure.

    model = structure[0]
    dssp = DSSP(model, in_file = pdb_file, dssp='dssp')

    if len(dssp.keys()) == 0:
        raise ValueError('DSSP could not be calculated for the given structure. Please check the input file.')
        return

    #First step: determining which residues in the training file are from a binding site
    #Iterating to find heteroatoms and storing them.

    if binding: 
        HETlist = []
        for model_instance in structure.get_iterator():
            for chain_instance in model_instance.get_iterator():
                for residue_instance in chain_instance.get_iterator():
                    if residue_instance.get_id()[0][0] == 'H':
                        HETlist.append(residue_instance)
        
        #Creating a dictionary to store the desired features.
                        
        dicc = {'protein': [], 'chainID': [], 'residueID': [], 'residue_instance': [], 'binding': [], 'residue': [], 'kd_hydrophobicity': [], 'rose_hydrophobicity': [],
                'charge': [], 'polarizability': [], 'molar_mass': [], 'volume_vdw': [], 'pI': [], 'steric_parameter': [], 'secondary_structure': [],
                'ASA': [], 'phi': [], 'psi': [], 'NH_O_1_energy': [], 'O_NH_1_energy': [], 'NH_O_2_energy': [], 'O_NH_2_energy': [], 'neighbor_density': [],
                'ASA_Zscore': []}
    else: 

        #Creating a dictionary to store the desired features.

        dicc = {'protein': [], 'chainID': [], 'residueID': [], 'residue_instance': [], 'residue': [], 'kd_hydrophobicity': [], 'rose_hydrophobicity': [],
                'charge': [], 'polarizability': [], 'molar_mass': [], 'volume_vdw': [], 'pI': [], 'steric_parameter': [], 'secondary_structure': [],
                'ASA': [], 'phi': [], 'psi': [], 'NH_O_1_energy': [], 'O_NH_1_energy': [], 'NH_O_2_energy': [], 'O_NH_2_energy': [], 'neighbor_density': [],
                'ASA_Zscore': []}
        
    #Iterating trough the structure and getting those residues that are not heteroatoms.
        
    for model_instance in structure.get_iterator():
        for chain_instance in model_instance.get_iterator():
            chainID = chain_instance.get_id()
            for residue_instance in chain_instance.get_iterator():
                residueID = residue_instance.get_id() 
                if residueID[0][0] == ' ':

                    # Check if the residue has structure information in the DSSP result; otherwise skip it

                    dssp_key = (chainID, residueID)
                    try:
                        dssp_info = dssp[dssp_key]
                    except KeyError: 
                        continue

                    #Storing the residue instance in the dictionary.

                    dicc['residue_instance'].append(residue_instance)

                    #If binding is set to true, determining if a residue is a binding site (has a minimum distance lower than  4 with at least one heteroatom)

                    if binding:
                        binding_res = 0
                        for HETres in HETlist:
                            dist = minimum_distance(HETres, residue_instance)
                            if dist <= 4:
                                binding_res = 1

                        #Appending 1 to the dictionary if a residue is a binding site and 0 if it is not.
                        
                        dicc['binding'].append(binding_res) 

                    #Getting the name and the residue  letter.
                        
                    residue = residue_instance.get_resname()
                    residue_letter = aa.aa_names.get(residue)

                    #Appending the pdb_code, the chain ID, the residue ID and the residue_letter to the dictionary.

                    dicc['protein'].append(pdb_code)
                    dicc['chainID'].append(chainID)
                    dicc['residueID'].append(residueID)
                    dicc['residue'].append(residue_letter)

                    # Appending the below features for our residue to the dictionary (they are obtained for the amino_acid_info module).

                    dicc['kd_hydrophobicity'].append(aa.kd_hydrophobicity.get(residue_letter)) #Kd_hydrophobicity.
                    dicc['rose_hydrophobicity'].append(aa.rose_hydrophobicity.get(residue_letter)) #rose_hydrophobicity.
                    dicc['charge'].append(aa.residue_charge.get(residue_letter)) #Charge
                    dicc['polarizability'].append(aa.polarizability.get(residue_letter)) #Polarizability.
                    dicc['molar_mass'].append(aa.weights.get(residue_letter)) #Molar mass.
                    dicc['volume_vdw'].append(aa.volume_van_der_waals.get(residue_letter)) #Volume van der Waals.
                    dicc['pI'].append(aa.pI.get(residue_letter)) #Isoelectric point.
                    dicc['steric_parameter'].append(aa.steric_parameter.get(residue_letter)) #Steric parameter

                    #Accessing the info of our residue stored in the dssp table and adding it to the dictionary. 
                    dicc['secondary_structure'].append(dssp_info[2]) #Secondary structure
                    dicc['ASA'].append(dssp_info[3]) #relative ASA
                    dicc['phi'].append(dssp_info[4]) #phi angle
                    dicc['psi'].append(dssp_info[5]) #psi angle
                    dicc['NH_O_1_energy'].append(dssp_info[7])  #NH–>O_1_energy
                    dicc['O_NH_1_energy'].append(dssp_info[9]) #O–>NH_1_energy
                    dicc['NH_O_2_energy'].append(dssp_info[11]) #NH–>O_2_energy
                    dicc['O_NH_2_energy'].append(dssp_info[13]) #O–>NH_2_energy

                    #Calculating the neighbor density and the ASA_Z and adding it to the dictionary.

                    density_info = calculate_density_and_convexity(structure, residue_instance, 10, dssp_info[3], dssp)
                    dicc['neighbor_density'].append(density_info[0]) #Neighbor density.
                    dicc['ASA_Zscore'].append(density_info[1]) #ASA_Zscore.
    
    #Converting the dictionary into a pandas dataframe. 
                    
    df = pd.DataFrame(dicc) 

    #Returning the dataframe. 
    
    return(df)