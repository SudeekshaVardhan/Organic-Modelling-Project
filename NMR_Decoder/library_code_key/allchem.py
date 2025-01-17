import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem

# Library provides additional functionalities

'''
Functionality:
- Molecular descriptors
- Molecular optimization
- Conformer generation
- Molecular reaction simulation
- Molecular mechanics
- Quantum mechanical calculations
'''

# read a SMILES string
molecule = Chem.MolFromSmiles("CCO")

# generate a set of low energy conformers
AllChem.EmbedMolecule(molecule)
