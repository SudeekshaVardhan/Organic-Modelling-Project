from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

class Input:
    def __init__(self):
        self.molecule = None

    def takeInput(self):
        self.molecule = input("Enter molecule: ")
    
    def returnInput(self):
        print(self.molecule)

    def readInput(self):
        reiterate = ""
        # identifies the functional groups
        for char in self.molecule:
            # check for benzene (how to do other cyclical molecules? learn how to simplify)
            if "H6" or "H5" in self.molecule:
                temp = self.molecule
                words = temp.split("C6")
                self.molecule = "c1ccccc1".join(words)
                                           
            if "COOH" in self.molecule:
                temp1 = self.molecule
                words1 = temp1.split("COOH")
                print(words1)
            else:
                break
        # remover step (could turn this into a method potentiallly)
        for char in self.molecule:
            # checks for additional H (not read in SMILES format)
            if char == 'H':
                continue
            # checks for numeric values (SMILES assumes numbers)
            elif char == str(2) or char == str(3) or char == str(5) or char == str(6):
                continue
            # checks for dashes that indicate single bonds in the molecule
            elif char == "-":
                continue
            else:
                reiterate += char

        print(reiterate)

        m = Chem.MolFromSmiles(reiterate)
        img = Draw.MolToImage(m)

        img.save('new.png')
        img.show()

            
            

