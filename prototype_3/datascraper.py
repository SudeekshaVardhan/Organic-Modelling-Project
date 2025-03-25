# Necessary Libraries
import vtk
import vtkmodules.vtkInteractionStyle
import vtkmodules.vtkRenderingOpenGL2

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pubchempy as pcp

from mp_api.client import MPRester 

import os


# This class is to search and process data from the pubchem database
class Datascraper:
    
    molecule = None
    numAtoms = []

    def __init__(self):
        pass

    def datascraping(self, molecule): 
        self.molecule = molecule
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
        # Returns the number of atoms in a molecule
        data = pcp.get_compounds(self.molecule, 'name')[0] # Checks for the molecule in the PubChem database
        mol = Chem.MolFromSmiles(data.canonical_smiles)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol) 
        print(mol.GetNumAtoms())

    def returnNames(self, mol_name):
        try:
            if not mol_name: # Checks if there is a valid molecule entry
                return "No molecule found. Try again."
            
            data = pcp.get_compounds(mol_name, "name")[0]
            compounds = list(set(data.synonyms))
            
            nicer = [s for s in compounds if not any (char.isdigit() for char in s[:5])]
            nicer = sorted(nicer, key=len)
    
            for x in nicer:
                print(x)
                print("")

        except Exception as e:
                print(e)
                return

    def getIsomers(self, mol_name):
        try:
            if not mol_name:
                return "No"
            data = pcp.get_compounds(mol_name, "name")[0]
            st1 = data.defined_atom_stereo_count
            st2 = data.defined_bond_stereo_count
            if st1 > 0 and st2 > 0:
                return True            
            for compound in data:
                print("o")
        except Exception as e:
            print(e)
            return

class CSVScrape:
    # Uses the Materials Project database - API Key required to get access (programmer - not sure if user needs it yet)
    # My account for API is connected to my GitHub account

    def __init__(self):
        self.network = None
        self.API_KEY = "CcWdfjf4GIvmnrTL4QkPOGPuF9JauIuo" # My unique API Key from the account
        
    def getID(self, name):
        with MPRester(self.API_KEY) as mpr:
            try:
                # Query the database using the material's formula
                structure = mpr.get_structure_by_formula(network)
                
                # Fetch the CIF data
                cif_data = structure.to(fmt="cif")
                
                # Write CIF to a file
                with open(f"{network}.cif", "w") as file:
                    file.write(cif_data)
                print(f"Downloaded CIF file for {network}!")
            except Exception as e:
                print(f"Error: {e}")

    def getNetwork(self, network):
        self.network = network
        mpid = self.getID(self.network)
        
        if mpid:
            try:
                with MPRester(self.API_KEY) as mpr:
                    # Fetch the CIF structure using the MP-ID
                    structure = mpr.get_structure_by_task_id(mpid)
                    
                    # Save the CIF structure to a file
                    cif = structure.to(fmt="cif")
                    with open("network.cif", "w") as file:
                        file.write(cif)
                    
                    print("Downloaded CIF file for 3D modeling!")
            except Exception as e:
                print(f"Failed to fetch CIF file: {e}")
        else:
            print("Failed to retrieve MP-ID, cannot fetch CIF file.")

class Modelling (Datascraper, CSVScrape):
    
    newMol = vtk.vtkMolecule()
    
    def __init__(self):
        super().__init__()
    
    def createMol(self):
        
        self.newMol.Initialize() # Prevents the molecules from stacking on top of each other
       
        try:
            if self.molecule == None:
                print("Error: No input for molecule")
            stuff = pcp.get_compounds(self.molecule, "name")[0]
            if not stuff:
                print("Error: Mol not found")
            mol = Chem.MolFromSmiles(stuff.canonical_smiles)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)

            for atom in mol.GetAtoms():
                self.newMol.AppendAtom(atom.GetAtomicNum(), mol.GetConformer().GetAtomPosition(atom.GetIdx()))

            for bond in mol.GetBonds():
                self.newMol.AppendBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())
        except Exception as e:
            print(e)
    
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

