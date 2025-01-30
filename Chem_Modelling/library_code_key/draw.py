from rdkit import Chem
from rdkit.Chem import Draw

# Pillow is needed to draw images
from PIL import Image

'''
SMILES syntax to know

Defaults
- C - Methane (CH4)
- O - Water (H2O)
- S - Hydrogen sulfide (H2S)
- N - Ammonia (NOT AMMONIUM - NH3)
- P - Phosphine (PH3)
- Cl - Hydrogen Chloride (HCl)

When making more complex molecules, H is usually not shown

Elements are shown in brackets: [Au] is elemental gold

Bonds:
- Single bonds usually not written (CC is CH3CH3 (ethane))
- Double bonds represented by = (C=O is CH2=O)
- Triple bonds represented by # ()
'''

# Draw the image
m = Chem.MolFromSmiles('CCO')
img = Draw.MolToImage(m)

# Save the image as a skeletal structure (2D rendering)
img.save("molecule.png")
img.show()