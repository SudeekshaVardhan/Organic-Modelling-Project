from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

class Input:
    
    def __init__(self, molecule):
        self.molecule = None
    
    def returnInput(self):
        print(self.molecule)

    def readInput(self):
        reiterated = ''
        nums = "1234567890"
        for char in self.molecule:
            # checks for additional H (not read in SMILES format)
            if char == 'H':
                continue
            # checks for numeric values (SMILES assumes numbers)
            elif char == str(2) or char == str(3):
                continue
            elif char == "-":
                continue
            else:
                reiterated += char
        print(reiterated)

        m = Chem.MolFromSmiles(reiterated)
        img = Draw.MolToImage(m)

        img.save('new.png')

            
            

