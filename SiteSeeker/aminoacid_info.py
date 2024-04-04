# Amino acid name conversion
aa_names = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# Kyte-Doolittle hydrophobicity scale
kd_hydrophobicity = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# Rose hydrophobicity scale: average area of buried amino acids in globular proteins
rose_hydrophobicity = {
    'A': 0.74, 'R': 0.64, 'N': 0.63	, 'D': 0.62, 'C': 0.91,
    'Q': 0.62, 'E': 0.62, 'G': 0.72, 'H': 0.78, 'I': 0.88,
    'L': 0.85, 'K': 0.52, 'M': 0.85, 'F': 0.88, 'P': 0.64,
    'S': 0.66, 'T': 0.70, 'W': 0.85, 'Y': 0.76, 'V': 0.86
}

# Amino acid charge at pH 7.4 (Histidine: Positive, 10%; Neutral, 90%)
residue_charge = {
    'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0,
    'Q': 0, 'E': -1, 'G': 0, 'H': 0.1, 'I': 0,
    'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0,
    'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0
}

# Polarizability
polarizability = {
    'A': 0.05, 'R': 0.29, 'N': 0.13, 'D': 0.11, 'C': 0.13,
    'Q': 0.18, 'E': 0.15, 'G': 0.00, 'H': 0.23, 'I': 0.19,
    'L': 0.19, 'K': 0.22, 'M': 0.22, 'F': 0.29, 'P': 0.00,
    'S': 0.06, 'T': 0.11, 'W': 0.41, 'Y': 0.30, 'V': 0.14
}

# Amino acids molar mass
weights = {
    'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1,
    'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 
    'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 
    'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2,
    'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19
}

# Volume normalized by Van der Waals volume
volume_van_der_waals = {
    'A': 1.00, 'G': 0.00, 'V': 3.00, 'L': 4.00,
    'I': 4.00, 'F': 5.89, 'Y': 6.47, 'W': 8.08,
    'T': 2.60, 'S': 1.60, 'R': 6.13, 'K': 4.77,
    'H': 4.66, 'D': 2.78, 'E': 3.78, 'N': 2.95,
    'Q': 3.95, 'M': 4.43, 'P': 2.72, 'C': 2.43
}

# isoelectric point
pI = {
    'A': 6.11, 'G': 6.07, 'V': 6.02, 'L': 6.04, 
    'I': 6.04, 'F': 5.67, 'Y': 5.66, 'W': 5.94, 
    'T': 5.60, 'S': 5.70, 'R': 10.74, 'K': 9.99, 
    'H': 7.69, 'D': 2.95, 'E': 3.09, 'N': 6.52, 
    'Q': 5.65, 'M': 5.71, 'P': 6.80, 'C': 6.35
}

# steric parameter: measure of the spatial requirements or bulkiness of the side chain of the amino acid residue
steric_parameter = {
    'A': 1.28, 'G': 0.00, 'V': 3.67, 'L': 2.59,
    'I': 4.19, 'F': 2.94, 'Y': 2.94, 'W': 3.21,
    'T': 3.03, 'S': 1.31, 'R': 2.34, 'K': 1.89,
    'H': 2.99, 'D': 1.60, 'E': 1.56, 'N': 1.60,
    'Q': 1.56, 'M': 2.35, 'P': 2.67, 'C': 1.77
}