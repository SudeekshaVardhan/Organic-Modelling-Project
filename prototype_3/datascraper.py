from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import numpy as np
import pubchempy as pcp
import plotly
import plotly.graph_objs as go

compound_list = [] # static variable, refreshes each time program is run. while program is open, resets


# This class is to search and process data from the pubchem database
class Datascraper:
    
    molecule = None
    numAtoms = []

    def __init__(self):
        pass

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
        embed = mol.GetConformer().GetPositions()

    # "Getter" methods
    # Returns data using methods from pcp library
    # Work in progress!!!
    '''
        What it should return:
        - Molecular formula
        - Properties
    '''
    def returnNumAtoms(self):
        data = pcp.get_compounds(self.molecule, 'name')[0]
        mol = Chem.MolFromSmiles(data.canonical_smiles)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        print(mol.GetNumAtoms())

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
                mol_name = pcp.get_compounds(mol, 'name')[0]
                data = pcp.get_properties(mol_name)
                print(data)
            else:
                print(pcp.NotFoundError)

    def returnMolWeight(self, mol_name):
        for mol in compound_list:
            if mol == mol_name:
                data = pcp.get_compounds(mol, 'name')[0]
                print(data.molecular_weight)
            else:
                print(pcp.NotFoundError)

class Modelling (Datascraper):
    def __init__(self):
        super().__init__()

    
