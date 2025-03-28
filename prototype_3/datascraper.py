# Necessary Libraries

# Chem Libraries
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pubchempy as pcp
from  thermo.chemical import Chemical # Used for thermometric data

# Modelling Libraries
import vtk
import vtkmodules.vtkInteractionStyle
import vtkmodules.vtkRenderingOpenGL2
import os


# This class is to search and process data from the pubchem database
class MolSearch:

    molecule = None # Molecule is accessible to both the parent and child classes (public variable that is accessible)

    def __init__(self):
        pass

    def datascraping(self, molecule): 

        self.molecule = molecule
        data = pcp.get_compounds(self.molecule,'name')[0] # Get the first matching entry from PubChem (if it exists)
        if not data: # If molecule does not exist in the PubChem database
            print(f"Error: {self.molecule} not found in PubChem")
            return
        
        # Get the molecular strucuture from the SMILES values stored in PubChem
        # PubChempy stories SMILES, and using this you can get the mol ()
        mol = Chem.MolFromSmiles(data.isomeric_smiles)
        mol = Chem.AddHs(mol)

        # If molecule is empty
        if mol is None:
            print(f"Error: Invalid molecule")
            return

        # Save the molecule to the same directory as the file path
        script_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(script_dir, "mol.png")
        img = Draw.MolToFile(mol, file_path, size=(200,200))

    '''
    "Getter" methods - 2 total
    Returns data using methods from pcp library
        What it should return:
        - Alternate Names
        - Properties
            - Log P (used for acids/bases)
            - Compound charge
            - IUPAC Name (great for practice for Structure 3)
            - Boiling Point
            - Melting Point (both these are used when looking at IMF)
            - Enthalpy (Energy stored in one mol of an atom. The more negative, the more energy)
            - Entropy (How "ordered" a molecule is - the more negative, the less chaotic)
            - Gibbs Free Energy (whether a reaction is spontaneous - used with other compounds)
    '''
    def returnNames(self, mol_name):
        try:
            if not mol_name: # Checks if there is a valid molecule entry
                return "No molecule found. Try again."
            
            data = pcp.get_compounds(mol_name, "name")[0] # Get compound from PubChem (first entry)
            compounds = list(set(data.synonyms))
            
            # This is the list of names returned - more user-friendly and ordered from shortest to longest
            nicer = [s for s in compounds if not any (char.isdigit() for char in s[:5])] # Filter out any number strings
            nicer = sorted(nicer, key=len)
    
            return nicer
        
        # Returns any uncaught errors in the terminal
        except Exception as e:
                print(e)
                return

    # Properties found in IB data booklet for select compounds (things that IB Chem students need)
    # For enthalpy, the temperature is fixed at 298 (the standard values in the databooklet for compounds)

    def getProperties(self, mol_name, param=None):
        try:
            if not mol_name: # Ensure that user is passing in a valid molecule
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

            if param: # If param exists, then get the info at the parameter
                return properties.get(param, f"Error: '{param}' is not a valid property")

            # If no specified parameter given, return the entire dictionary
            return properties
        
        except Exception as e:
            return f"Error: {e}"  # Return error message
        
class MolModelling(MolSearch):
    
    newMol = vtk.vtkMolecule() # Also accessible outside of the class - public variable
    
    def __init__(self):
        # The MolModelling is a child of MolSearch. This is so that the 3D rendering capabilities can 
        super().__init__()
    
    def createMol(self):
        
        self.newMol.Initialize() # Prevents the molecules from stacking on top of each other
       
        try:
            if self.molecule == None:
                print("Error: No input for molecule")
            stuff = pcp.get_compounds(self.molecule, "name")[0] # same as the datascraper in the previous class
            if not stuff:
                print("Error: Mol not found")
            mol = Chem.MolFromSmiles(stuff.isomeric_smiles)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)


            # Add atoms based on the atomic number and (x,y,z) coordinates
            # The atomic num also determines the color of the molecule
            for atom in mol.GetAtoms():
                self.newMol.AppendAtom(atom.GetAtomicNum(), mol.GetConformer().GetAtomPosition(atom.GetIdx()))

            # For each bond found, add it to the VTK molecule using the beginning, end, and type of bond between atoms
            for bond in mol.GetBonds():
                self.newMol.AppendBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())
        except Exception as e:
            print(e)

    def renWin(self):

        # Objects used to map the molecule in 3D
        newMapper = vtk.vtkMoleculeMapper()
        newActor = vtk.vtkActor()
        newRenderer = vtk.vtkRenderer()
        newRenWin = vtk.vtkRenderWindow()

        # Assign values to the newMol object
        self.createMol()

        newMapper.SetInputData(self.newMol)
        newActor.SetMapper(newMapper)
        newRenderer.AddActor(newActor)
        newRenWin.AddRenderer(newRenderer)
        newRenWin.SetSize(500,800)

        # Allows for user interaction with the molecule (and it holds it open until user closes out of window)
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(newRenWin)

        newRenWin.Render()
        interactor.Initialize()
        interactor.Start()
