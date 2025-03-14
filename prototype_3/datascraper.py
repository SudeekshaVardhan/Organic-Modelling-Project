# Necessary Libraries
import vtk
import vtkmodules.vtkInteractionStyle
import vtkmodules.vtkRenderingOpenGL2

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pubchempy as pcp

import os

compound_list = [] # static variable, refreshes each time program is run. while program is open, resets

# This class is to search and process data from the pubchem database
class Datascraper:
    
    molecule = None
    numAtoms = []

    def __init__(self):
        pass

    def datascraping(self, molecule): 
        self.molecule = molecule
        compound_list.append(self.molecule)
        data = pcp.get_compounds(self.molecule,'name')[0]     
        if not data:
            print(f"Error: Molecule not found in PubChem")
            return
        
        mol = Chem.MolFromSmiles(data.canonical_smiles)
        mol = Chem.AddHs(mol)

        if mol is None:
            print(f"Error: Invalid molecule")
            return
        
        # Save in the same directory as program
        script_dir = os.path.dirname(os.path.abspath(__file__))  
        file_path = os.path.join(script_dir, "mol.png")  # Save inside the same directory

        # Save image
        img = Draw.MolToFile(mol, file_path, size=(200,200))
        
        # Check it was saved
        print(f"mol.png successfully saved at: {file_path}")

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
    # Work in progress - redo so not dependent on the list object
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
    # Used to read networks from other databases. Should access both pubchem and a seperate database or file for coordinates.
    # Incomplete - additional work must be done on this
    def __init__(self):
        pass


class Modelling (Datascraper):
    
    newMol = vtk.vtkMolecule()
    
    def __init__(self):
        super().__init__()
    
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
    
    def createNet(self):
        pass

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