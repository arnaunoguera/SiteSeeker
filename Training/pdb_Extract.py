from rcsbsearchapi.const import CHEMICAL_ATTRIBUTE_SEARCH_SERVICE, STRUCTURE_ATTRIBUTE_SEARCH_SERVICE
from rcsbsearchapi.search import AttributeQuery
import urllib.request
from pathlib import Path

with open('ligands.txt', 'r') as file:
    for line in file:

        chemID = line.strip()
        query = AttributeQuery("rcsb_chem_comp_container_identifiers.comp_id", "exact_match", chemID, CHEMICAL_ATTRIBUTE_SEARCH_SERVICE)
        # Call the query to execute it, keep only the first protein ID (to avoid overfitting to a certain type of compound)
        for assemblyid in query("assembly"):
            pdb_id = assemblyid.split('-')[0]
            urllib.request.urlretrieve('http://files.rcsb.org/download/'+pdb_id+'.pdb', 'pdb_files/' + pdb_id + '.pdb')
            print("Successfully retrieved " + pdb_id + " for compound " + chemID)
            break
