from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pubchempy as pcp

from PIL import Image

# Fetch the compound from PubChem
compound_list = pcp.get_compounds('Glucose','name')

# Access the first compound in the list
compound = compound_list[0]
print(compound)

canonical_smiles = compound.canonical_smiles
img = Draw.MolToImage(str(canonical_smiles))

img.save("hello.png")
