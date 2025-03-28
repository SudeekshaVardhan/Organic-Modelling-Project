
# 2D rendering libraries + PubChempy for additional info:
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pubchempy as pcp
from thermo.chemical import Chemical

# 3D rendering libraries (get and read CIF files and other data)
import requests
import json
import os
import vtk
import numpy as np
import time
from ase import io
from mp_api.client import MPRester
from pymatgen.io.cif import CifWriter


class NetworkSearch:
    # Uses both the COD database and the MP database
    # COD is used for getting the formula of the compound entered (as the user is prompted to enter the name)
    # MP is used to get the CIF files of each molecule
    # My account for API is connected to my GitHub account

    def __init__(self):
        self.network = None # Takes the user input and resets it
        self.formula = None  # Gets it from get formula method (uses COD to get formula). Input this into MP database
        self.struct = None
        
    def getCIF(self, name):
        self.network = name
        # Searches the database by compound name (Crystallography Open Database) to get the formula.
        url = f"https://www.crystallography.net/cod/result?text={self.network}&format=json" 
        response = requests.get(url)

        if response.status_code == 200: # Check if database can be accessed
            try:
                data = response.json() # Get list of json that have the same entry
                if isinstance(data, list) and data:
                    form = data[0].get("formula", "") # Get the first entry in the database that matches entered formula
                    self.formula = form.replace(" ", "").strip("-")
                    if not self.formula:
                        return "Formula doesn't exist"
                else:
                    return "Error: List does not exist"
                
            except json.JSONDecodeError:
                return f"Failed to read JSON file: {json.JSONDecodeError}"
        
        else: # Reason why it doesn't read database (Database error)
            return f"Error: Database is currently not working (Error Code {response.status_code})"
        
      
        print(self.formula)
      
        # Get CIF File from MP Rester

        myAPI = "EMhu1YM4AQ883JBFbW7wl19zyI0NkqTA"
        
        if not self.formula:
            return "Formula doesn't exist. Try again"
        
        save = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(save, "networkCoord.cif")
        
        try:
            with MPRester(api_key=myAPI) as rester:

                results = rester.materials.summary.search(
                     formula=self.formula,
                     fields=["material_id"]    
                )    

            id = results[0].material_id # First result 
            if id:
                print(f"{id}")
            else:
                print(f"{id} doesn't exist")
            self.struct = rester.get_structure_by_material_id(id) # Get structure from ID (Unique to MP database - cant use COD ID)

            writer = CifWriter(self.struct)
            writer.write_file(file_path)
            time.sleep(1) # Account for lag in updating network coords

            if os.path.exists(file_path):
                return f"File saved successfully at {file_path}!"
            else:
                return "Oh no!"

        except Exception as e:
            return f"Error: {e}" # Return any errors while downloading the data from file
        
    # Fixed values for temperature @ 298 K (similar method in the MolSearch class)
    def getProperties(self, mol_name, param=None):
        try:
            if not mol_name:
                return "No molecule entered"
            
            # Get data from PubChem
            data = pcp.get_compounds(mol_name, "name")[0]
            
            # Object of the thermo.chemical library. Returns thermodynamic data on select compounds.
            chem = Chemical(mol_name)

            # Create a dictionary for each
            properties = {
                "logP": getattr(data, "xlogp", "N/A"),
                "charge": getattr(data, "charge", "N/A"),
                "iupac": getattr(data, "iupac_name", "N/A"),
                "boiling_point": getattr(data, "boiling_point", "N/A"),
                "melting_point": getattr(data, "melting_point", "N/A"),
                "specific_heat_capacity": getattr(chem, "Cp", "N/A"),
                "enthalpy": getattr(chem, "Hf", "N/A"),
                "entropy": getattr(chem, "Sf", "N/A"),
                "gibbs_free_energy": getattr(chem, "Gf", "N/A")
            }
        
            if param: # If param exists, then get the info at the 
                return properties.get(param, f"Error: '{param}' is not a valid property")

            return properties
        except Exception as e:
            return str(e)
    
    def getNetImg(self, name):
        # Gets info from the PubChempy database instead of COD or MP Rester to generate 2D image. 
        # Properties and other names generated using this database (like uses, etc.)
        self.network = name
        data = pcp.get_compounds(self.network,'name')[0]     
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
        file_path = os.path.join(script_dir, "network.png")  # Save inside the same directory

        # Save image
        img = Draw.MolToFile(mol, file_path, size=(200,200))
        
        # Check it was saved
        print(f"mol.png successfully saved at: {file_path}")

class NetModel(NetworkSearch):

    poly = vtk.vtkPolyData()

    def __init__(self):
        super().__init__()
        self.point = vtk.vtkPoints()
        self.bonds = vtk.vtkCellArray()
        self.types = []
    
    def newModel (self):
        # Similar to the VTK of the molecule
        # To prevent stacking

        # Get the location of the model (where file is saved)
        path = os.path.join(os.path.dirname(__file__), "networkCoord.cif")

        if not os.path.exists(path):
            return FileNotFoundError(f"No file at {path}")        
        
        atoms = io.read(path)

        positions = atoms.positions
        types = atoms.get_chemical_symbols()
        print(f"types {types}")
        print(f"positions: {positions}")

        for i, position in enumerate(positions):
            self.point.InsertNextPoint(position[0], position[1], position[2])
            self.types.append(types[i])
            print(f"Atom added at {position}")

        thres = 2.0
        num = len(positions)

        for i in range(num):
            for j in range(i+1, num):
                dist = np.linalg.norm(positions[i] - positions[j])
                if dist < thres:
                    self.bonds.InsertNextCell(2)
                    self.bonds.InsertCellPoint(i)
                    self.bonds.InsertCellPoint(j)
                    
    def netRenWin(self): 

        self.poly.Initialize()
        
        self.newModel()
        self.poly.SetPoints(self.point)
        self.poly.SetLines(self.bonds)

        newMapper = vtk.vtkPolyDataMapper()
        newMapper.SetInputData(self.poly)

        newActor = vtk.vtkActor()
        newActor.SetMapper(newMapper)

        newRenderer = vtk.vtkRenderer()
        newRenderer.AddActor(newActor)

        newRenWin = vtk.vtkRenderWindow()
        newRenWin.AddRenderer(newRenderer)

        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(newRenWin)

        newRenWin.Render()
        interactor.Start()
