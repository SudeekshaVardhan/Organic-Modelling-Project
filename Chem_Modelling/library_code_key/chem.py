from rdkit import Chem

# Basic chem informatic functionality for working with molecular structures
# Reads molecule and handles different file formats

'''
Functionality
- Reading and writing mols in various file formats
- Calculating properties
- Performing substructure searches
- Molecular Visualization
'''

mol = Chem.MolFromSmiles("CC(C(=O)O)N")
smiles = Chem.MolToSmiles(mol)

print(smiles)

