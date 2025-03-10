# Necessary Libraries
import vtk
import vtkmodules.vtkInteractionStyle
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderWindow, vtkRenderer

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
        if self.molecule.lower() == "graphene":
            stuff = "C12=C3C4=C5C6=C7C=8C(*)=C5C(=C(*)C4=C(C1=C(*)C9=C%10C2=C%11C%12=C3C6=C%13C%14=C7C=%15C=%16C%17=C%14C%18=C%19C%20=C%21C(=C%22C(C(*)=C%21C(*)=C(C%20=C(C%18=C(*)C%17=C(*)C%23=C(C(=C(C(C%16%23)=C(C%15C8*)*)*)*)*)*)*)=C(C=%24C(=C%11%22)C%10=C%25C(=C9*)C(*)=C(C(=C%25C%24*)*)*)*)C%12=C%19%13)*)*"
            mol = Chem.MolFromSmiles(stuff)
        elif self.molecule.lower() == "diamond":
            pass
        else:
            data = pcp.get_compounds(self.molecule,'name')[0]
            mol = Chem.MolFromSmiles(data.canonical_smiles)
        mol = Chem.AddHs(mol)
        img = Draw.MolToImage(mol)
        img.show() # Display image in local image editor (remove in later versions, or push to seperate GUI)

    def returnSMILES(self):
        data = pcp.get_compounds(self.molecule,'name')[0]
        mol = Chem.MolFromSmiles(data.canonical_smiles)
        mol = Chem.AddHs(mol)
        print(mol)

    def returnCoords(self):
        data = pcp.get_compounds(self.molecule, 'name')[0]
        if not data:
            return pcp.NotFoundError()
        mol = Chem.MolFromSmiles(data.canonical_smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        embed = mol.GetConformer().GetPositions()
        print(embed)

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

class CSVScrape:
    def __init__(self):
        pass


class Modelling (Datascraper):
    
    newMol = vtk.vtkMolecule()
    
    def __init__(self, data):
        super().__init__()
        self.data = Datascraper
    
    def createMol(self):
        stuff = pcp.get_compounds(self.molecule, "name")[0]
        mol = Chem.MolFromSmiles(stuff.canonical_smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)

        for atom in mol.GetAtoms():
            self.newMol.AppendAtom(atom.GetAtomicNum(), mol.GetConformer().GetAtomPosition(atom.GetIdx()))

        for bond in mol.GetBonds():
            self.newMol.AppendBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())
    


    def renWin(self):

        newMapper = vtk.vtkMoleculeMapper()
        newActor = vtk.vtkActor()
        newRenderer = vtk.vtkRenderer()
        newRenWin = vtk.vtkRenderWindow()
        self.createMol()

        newMapper.SetInputData(self.newMol)
        newActor.SetMapper(newMapper)
        newRenderer.AddActor(newActor)
        newRenWin.AddRenderer(newRenderer)
        newRenWin.SetSize(500,800)

        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(newRenWin)

        newRenWin.Render()
        interactor.Initialize()
        interactor.Start()