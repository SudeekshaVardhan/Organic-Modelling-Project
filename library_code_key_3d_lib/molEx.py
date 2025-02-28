# Necessary Libraries
import vtk
import vtkmodules.vtkInteractionStyle
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderWindow, vtkRenderer

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

import numpy as np
import pubchempy as pcp
import plotly
import plotly.graph_objs as go

info = input("Enter mol:")
data = pcp.get_compounds(info, "name")[0]
mol = Chem.MolFromSmiles(data.canonical_smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)
AllChem.UFFOptimizeMolecule(mol)
img = Draw.MolToImage(mol)
img.show()

newMol = vtk.vtkMolecule()

# Create points for the atom
for atom in mol.GetAtoms():
    newMol.AppendAtom(atom.GetAtomicNum(), mol.GetConformer().GetAtomPosition(atom.GetIdx()))

# Bond points
for bond in mol.GetBonds():
    newMol.AppendBond(bond.GetBeginAtomIdx(),bond.GetEndAtomIdx(), bond.GetBondType())

# Create points to map the molecule and generate an actor to render graphics
# Mapper
mapMol = vtk.vtkMoleculeMapper()
mapMol.SetInputData(newMol)

#Actor
molActor = vtkActor()
molActor.SetMapper(mapMol)

# Create the renderer
renderer = vtkRenderer()
renderer.AddActor(molActor)
renderer.SetBackground(0,0,0)

# Create the window
renWin = vtkRenderWindow()
renWin.AddRenderer(renderer)
renWin.SetSize(300,300)
print(renWin.GetRenderingBackend()) # Version used to run VTK (Open GL2)

# Create the interactor - allows user to rotate object in the window
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renWin)

renWin.Render()
interactor.Start()

'''
KEY:
- Carbon = gray
- Hydrogen = white
- Oxygen = red

'''
'''
Additional Notes:
- Check if it is possible to display the key on the vtk window
- Add functionality to tkinter (make this more pretty and user-friendly)
'''


