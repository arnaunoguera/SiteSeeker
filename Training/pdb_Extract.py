#Importing the necessary modules

from rcsbsearchapi.const import CHEMICAL_ATTRIBUTE_SEARCH_SERVICE, STRUCTURE_ATTRIBUTE_SEARCH_SERVICE
from rcsbsearchapi.search import AttributeQuery
import urllib.request
from pathlib import Path

#Opening a file with a list of ligands and iterating trough each ligand of the file.

with open('ligands.txt', 'r') as file:
    for line in file:
        chemID = line.strip()

        #Creating a query to search for the ligand ID in the RCSB database.

        query = AttributeQuery("rcsb_chem_comp_container_identifiers.comp_id", "exact_match", chemID, CHEMICAL_ATTRIBUTE_SEARCH_SERVICE)
        for assemblyid in query("assembly"):

            #For each ligand, storing only the ID of the first matching protein(this is done to avoid overfitting for certain ligands)

            pdb_id = assemblyid.split('-')[0]

            #Retrieving the pdb file for the stored protein from the RCSB website.
            
            urllib.request.urlretrieve('http://files.rcsb.org/download/'+pdb_id+'.pdb', 'pdb_files/' + pdb_id + '.pdb')
            print("Successfully retrieved " + pdb_id + " for compound " + chemID)
            break
