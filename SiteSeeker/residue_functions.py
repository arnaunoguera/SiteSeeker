from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBList
from Bio.PDB.DSSP import DSSP
import numpy as np
import pandas as pd
import SiteSeeker.aminoacid_info as aa

def get_distance(atom1, atom2):
    coord1 = np.array(atom1.get_coord())
    coord2 = np.array(atom2.get_coord())
    dist = np.sqrt(np.sum((coord1-coord2)**2))
    return dist

def minimum_distance(residue1, residue2):
    first = True
    for res1atom in residue1.get_iterator():
        for res2atom in residue2.get_iterator():
            dist = get_distance(res1atom,res2atom)
            if first == True or dist < min_dist:
                first = False
                min_dist = dist
    return min_dist

def calculate_density_and_convexity(structure, reference_residue, search_radius, ASA, dssp):
    neighbor_count = 0
    neighbor_list_ASA = []
    for residue2 in structure.get_residues():
        residue2ID = residue2.get_id() 
        if residue2ID[0][0] == ' ':
            distance = minimum_distance(reference_residue, residue2)
            if distance <= search_radius:
                neighbor_count += 1
                chainID = residue2.get_parent().get_id()
                dssp_key = (chainID, residue2ID)
                residue2ASA = dssp[dssp_key][3]
                neighbor_list_ASA.append(residue2ASA)
    # Calculate volume of a sphere with radius 'search_radius'
    sphere_volume = (4/3) * np.pi * (search_radius ** 3)
    
    # Calculate neighbor density
    neighbor_density = neighbor_count / sphere_volume

    #Calculate ASA_Z

    ASA_Z = ((ASA-np.mean(neighbor_list_ASA))/np.std(neighbor_list_ASA))
    
    return(neighbor_density, ASA_Z)

def residue_features(pdb_code, pdb_file, binding=True):
    p = PDBParser(PERMISSIVE=1, QUIET=1)
    
    structure = p.get_structure(pdb_code, pdb_file)

    model = structure[0]
    dssp = DSSP(model, in_file = pdb_file, dssp='dssp')

    #First step: determining which residues in the training file are from a binding site
    #Iterate to find heteroatoms and store them
    if binding: 
        HETlist = []
        for model_instance in structure.get_iterator():
            for chain_instance in model_instance.get_iterator():
                for residue_instance in chain_instance.get_iterator():
                    if residue_instance.get_id()[0][0] == 'H':
                        HETlist.append(residue_instance)
        dicc = {'protein': [], 'chainID': [], 'residueID': [], 'residue_instance': [], 'binding': [], 'residue': [], 'kd_hydrophobicity': [], 'rose_hydrophobicity': [],
                'charge': [], 'polarizability': [], 'molar_mass': [], 'volume_vdw': [], 'pI': [], 'steric_parameter': [], 'secondary_structure': [],
                'ASA': [], 'phi': [], 'psi': [], 'NH_O_1_energy': [], 'O_NH_1_energy': [], 'NH_O_2_energy': [], 'O_NH_2_energy': [], 'neighbor_density': [],
                'ASA_Zscore': []}
    else: 
        dicc = {'protein': [], 'chainID': [], 'residueID': [], 'residue_instance': [], 'residue': [], 'kd_hydrophobicity': [], 'rose_hydrophobicity': [],
                'charge': [], 'polarizability': [], 'molar_mass': [], 'volume_vdw': [], 'pI': [], 'steric_parameter': [], 'secondary_structure': [],
                'ASA': [], 'phi': [], 'psi': [], 'NH_O_1_energy': [], 'O_NH_1_energy': [], 'NH_O_2_energy': [], 'O_NH_2_energy': [], 'neighbor_density': [],
                'ASA_Zscore': []}
    for model_instance in structure.get_iterator():
        for chain_instance in model_instance.get_iterator():
            chainID = chain_instance.get_id()
            for residue_instance in chain_instance.get_iterator():
                residueID = residue_instance.get_id() 
                if residueID[0][0] == ' ':
                    dicc['residue_instance'].append(residue_instance)
                    if binding:
                        binding_res = 0
                        for HETres in HETlist:
                            dist = minimum_distance(HETres, residue_instance)
                            if dist <= 4:
                                binding_res = 1
                        dicc['binding'].append(binding_res) #append 1 if it's a binding site residue, 0 if not 
                    residue = residue_instance.get_resname()
                    residue_letter = aa.aa_names.get(residue)
                    dicc['protein'].append(pdb_code)
                    dicc['chainID'].append(chainID)
                    dicc['residueID'].append(residueID)
                    #print(pdb_code, chainID, residueID, sep = '\t')
                    dicc['residue'].append(residue_letter)
                    dicc['kd_hydrophobicity'].append(aa.kd_hydrophobicity.get(residue_letter))
                    dicc['rose_hydrophobicity'].append(aa.rose_hydrophobicity.get(residue_letter))
                    dicc['charge'].append(aa.residue_charge.get(residue_letter))
                    dicc['polarizability'].append(aa.polarizability.get(residue_letter))
                    dicc['molar_mass'].append(aa.weights.get(residue_letter))
                    dicc['volume_vdw'].append(aa.volume_van_der_waals.get(residue_letter))
                    dicc['pI'].append(aa.pI.get(residue_letter))
                    dicc['steric_parameter'].append(aa.steric_parameter.get(residue_letter))
                    #DSSP
                    dssp_key = (chainID, residueID)
                    dssp_info = dssp[dssp_key]
                    dicc['secondary_structure'].append(dssp_info[2]) #Secondary structure
                    dicc['ASA'].append(dssp_info[3]) #relative ASA
                    dicc['phi'].append(dssp_info[4]) #phi angle
                    dicc['psi'].append(dssp_info[5]) #psi angle
                    dicc['NH_O_1_energy'].append(dssp_info[7])  #NH–>O_1_energy
                    dicc['O_NH_1_energy'].append(dssp_info[9]) #O–>NH_1_energy
                    dicc['NH_O_2_energy'].append(dssp_info[11]) #NH–>O_2_energy
                    dicc['O_NH_2_energy'].append(dssp_info[13]) #O–>NH_2_energy
                    #NEIGHBOR DENSITY
                    density_info = calculate_density_and_convexity(structure, residue_instance, 10, dssp_info[3], dssp)
                    dicc['neighbor_density'].append(density_info[0])
                    dicc['ASA_Zscore'].append(density_info[1])
    df = pd.DataFrame(dicc)
    return(df)
