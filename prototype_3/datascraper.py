# Necessary Libraries

# Chem Libraries
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pubchempy as pcp
from  thermo.chemical import Chemical

# Modelling Libraries
import vtk
import vtkmodules.vtkInteractionStyle
import vtkmodules.vtkRenderingOpenGL2
import os


# This class is to search and process data from the pubchem database
class MolSearch:
    
    molecule = None
    numAtoms = []

    def __init__(self):
        pass

    def datascraping(self, molecule): 
        print(f"Scraping data for molecule: {molecule}")  # Debugging molecule name
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

        script_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(script_dir, "mol.png")
        
        print(f"Saving image to {file_path}")  # Debugging path
        img = Draw.MolToFile(mol, file_path, size=(200,200))
        print(f"mol.png successfully saved at: {file_path}")

    '''
    "Getter" methods
    Returns data using methods from pcp library

        What it should return:
        - Alternate Names
        - Properties
            - 

    '''
    def returnNames(self, mol_name):
        try:
            if not mol_name: # Checks if there is a valid molecule entry
                return "No molecule found. Try again."
            
            data = pcp.get_compounds(mol_name, "name")[0]
            compounds = list(set(data.synonyms))
            
            nicer = [s for s in compounds if not any (char.isdigit() for char in s[:5])]
            nicer = sorted(nicer, key=len)
    
            return nicer

        except Exception as e:
                print(e)
                return

    # Properties found in IB data booklet for select compounds (things that IB Chem students need)
    # For enthalpy, the temperature is fixed at 298 (the standard values in the databooklet for compounds)

    def getProperties(self, mol_name, param=None):
        try:
            if not mol_name:
                return "No molecule entered"
            
            # Get data from PubChem
            data = pcp.get_compounds(mol_name, "name")[0]
            
            # Object of the thermo.chemical library. Returns thermodynamic data on select compounds.
            chem = Chemical(str(mol_name))

            # Create a dictionary for each
            properties = {
                "logP": getattr(data, "xlogp", "N/A"),
                "charge": getattr(data, "charge", "N/A"),
                "iupac": getattr(data, "iupac_name", "N/A"),
                "boiling_point": getattr(data, "boiling_point", "N/A"),
                "melting_point": getattr(data, "melting_point", "N/A"),
                "specific_heat_capacity": chem.Cp,
                "enthalpy": chem.Hf,
                "entropy": chem.Sf,
                "gibbs_free_energy": chem.Gf
            }

            # Thermo data
        
            if param: # If param exists, then get the info at the 
                return properties.get(param, f"Error: '{param}' is not a valid property")

            return properties
        
        except Exception as e:
            return f"Error: {e}"  # Return error message
        
class MolModelling(MolSearch):
    
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
            mol = Chem.MolFromSmiles(stuff.isomeric_smiles)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)

            for atom in mol.GetAtoms():
                self.newMol.AppendAtom(atom.GetAtomicNum(), mol.GetConformer().GetAtomPosition(atom.GetIdx()))

            for bond in mol.GetBonds():
                self.newMol.AppendBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())
        except Exception as e:
            print(e)

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
