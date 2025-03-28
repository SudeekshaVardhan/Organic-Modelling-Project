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
        self.molecule = molecule
        data = pcp.get_compounds(self.molecule,'name')[0]     
        if not data:
            print(f"Error: Molecule not found in PubChem")
            return
        
        mol = Chem.MolFromSmiles(data.isomeric_smiles)
        mol = Chem.AddHs(mol)

        if mol is None:
            print(f"Error: Invalid molecule")
            return
        
        # Save in the same directory as program
        script_dir = os.getcwd()  
        file_path = os.path.join(script_dir, "mol.png")  # Save inside the same directory

        if os.path.exists(file_path):
            os.remove(file_path)

        # Save image
        img = Draw.MolToFile(mol, file_path, size=(200,200))

    def returnCoords(self):
        data = pcp.get_compounds(self.molecule, 'name')[0]
        if not data:
            return pcp.NotFoundError()
        mol = Chem.MolFromSmiles(data.isomeric_smiles)
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
        - Names
        - Properties
        - Isomers - Number of chiral centers, how many isomers

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
    
            return nicer

        except Exception as e:
                print(e)
                return

    # Properties found in IB data booklet for select compounds (things that IB Chem students need)
    # For enthalpy, the temperature is fixed at 298 (the standard values in the databooklet for compounds)
    def getProperties(self, mol_name, param):
        try:
            returnVal = param
            if not mol_name:
                return "No molecule entered"
            data = pcp.get_compounds(mol_name, "name")[0]
            
            logP = data.xlogp # Returns P value
            charge = data.charge # Get atom charge (if any)
            iupac = data.iupac_name # Get Atom Name
            boil = getattr(data, 'boiling_point', 'N/A')
            melt = getattr(data, 'melting_point', 'N/A')

            # Thermo data (get from the thermo.chemical library)
            chem = Chemical(mol_name)
            specHeatCap = chem.Cp(298)
            enthalpy = chem.H(298)
            entropy = chem.S(298)
            gibbs = chem.G(298)

            vals = [logP, charge, iupac, boil, melt, specHeatCap, enthalpy, entropy, gibbs]

            for val in vals:
                if param == val:
                    print(val)
                    return
                elif param == val and val == None:
                    print("N/A")
                    return
                else:
                    print("Request not found.")
                    return

        except Exception as e:
            return e
        
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

    