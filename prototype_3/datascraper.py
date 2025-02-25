from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import numpy as np
import pubchempy as pcp


compound_list = [] # static variable

import plotly
import plotly.graph_objs as go
# Stores the user's compounds (this is a static variable, so the user 
# can access any previous compound if they wish)

# This class is to search and process data from the pubchem database
class Datascraper:

    def __init__(self):
        self.molecule = None 

    def takeInput(self):
        user_molecule = input("Enter molecule name: ")
        self.molecule = user_molecule.lower()
        compound_list.append(self.molecule)
    
    def datascraping(self):
        data = pcp.get_compounds(self.molecule,'name')[0]
        mol = Chem.MolFromSmiles(data.canonical_smiles)
        mol = Chem.AddHs(mol)
        img = Draw.MolToImage(mol)
        img.show() # Display image in local image editor (remove in later versions, or push to seperate GUI)

    def getCoords(self):
        data = pcp.get_compounds(self.molecule, 'name')[0]
        mol = Chem.MolFromSmiles(data.canonical_smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        print(mol.GetNumAtoms())
        embed = mol.GetConformer().GetPositions()
        print(embed)

    # "Getter" methods
    # Returns data using methods from pcp library
    # Work in progress!!!
    '''
        What it should return:
        - Full array of compounds (print)
        - Molecular formula
        - Properties
    '''

    def returnData(self):
        for mol in compound_list:
            print(mol)

    def returnMF(self, mol_name):
        for mol in compound_list:
            if mol == mol_name:
                data = pcp.get_properties(mol)
                print(data)
            else:
                print(pcp.NotFoundError)

    def returnProperties(self, mol_name):
        for mol in compound_list:
            if mol == mol_name:
                ind = compound_list.index(mol)
                mol_name = pcp.get_compounds[0]
                mol = mol_name.get_canonicalsmiles
                data = pcp.get_properties(mol)
                print(data)
            else:
                print(pcp.NotFoundError)



