from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pubchempy as pcp

from PIL import Image

# Fetch the compound from PubChem
compound_list = pcp.get_compounds('Ethanoic Acid','name')

# Access the first compound in the list
compound = compound_list[0]
print(compound)

# Canonical smiles requires a list to work
canonical_smiles = compound.canonical_smiles

# Must convert from the SMILES string to mol object for image
mol = Chem.MolFromSmiles(canonical_smiles)
img = Draw.MolToImage(mol)

# Save and display image
img.save("hello.png")
img.show()