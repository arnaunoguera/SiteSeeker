# SiteSeeker

SiteSeeker is a tool for predicting binding sites in a protein structure based on a Support Vector Machine (SVM) model.
The SVM model considers a set of structural and physicochemical features and is trained on a dataset of protein structures with known binding sites.

For further theoretical explanations, a detailed tutorial and the analysis of examples, check out SiteSeeker_Report.pdf in the SiteSeeker folder. 

## Features

- Predicts binding sites in protein structures from a PDB file.
- Uses a SVM model trained on structural and physicochemical features.
- Outputs a PDB file containing the predicted binding residues.
- Generates a text file with information about the binding sites and the residues involved in them.
- Calculates accuracy, precision, and recall if the interaction with the ligand is already known.

## Usage
usage: 
```bash
SiteSeeker [-h] --input INPUT [--output OUTPUT] [--verbose] [--maxdist MAXDIST] [--knownbind]
```

Arguments:

  -h, --help            show the help message and exit

  --input INPUT, -i INPUT
                        Input file in PBD format (with extension .pdb)

  --output OUTPUT, -o OUTPUT
                        Desired path and base filename of the output files (without the file extension)

  --verbose, -v         Prints additional information during the execution

  --maxdist MAXDIST, -d MAXDIST
                        Maximum distance between the closest residues in two binding sites for these to be considered
                        separate binding sites (default: 4 Å)

  --knownbind, -k       The input PDB file contains the interaction with the ligand (with the flag HETATM) and I would
                        like to know the accuracy of SiteSeeker's prediction

## Output

- **PDB File**: A modified version of the input PDB file containing only the binding residues. This should be opened in Chimera alongside the original pdb file to visualize the results.
- **Text File**: Contains information about the predicted binding sites and the residues involved.

## Installation

pip install SiteSeeker-1.0.tar.gz

All required Python libraries should automatically be installed and SiteSeeker is added to the $PATH, so the program can be called in this way:
SiteSeeker [-h] --input INPUT [--output OUTPUT] [--verbose] [--maxdist MAXDIST] [--knownbind]

However, you should make sure DSSP is installed in your computer beforehand. The following code can be used to install DSSP on Ubuntu and make it available to SiteSeeker (obtained from https://ssbio.readthedocs.io/).

```bash
sudo apt-get update
sudo apt-get install dssp
sudo ln -s /usr/bin/mkdssp /usr/bin/dssp
```

## Examples

To predict binding sites in a protein structure:

```bash
SiteSeeker -i pdb_files/1HWK.pdb -o results/1HWK -v
```

To predict binding sites in a protein with a known ligand interaction and calculate accuracy, precision, and recall:

```bash
SiteSeeker -i pdb_files/1HWK.pdb -o results/1HWK -vk
```

## Support

You can contact us at arnau.noguera01@estudiant.upf.edu and aina.vaquer01@estudiant.upf.edu 