
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Fingerprints

# example molecule in SMILES format

smiles = "CCO"
mol = Chem.MolFromSmiles(smiles)

# generate the Morgan fingerprint (radius = 2, 2048 bits)