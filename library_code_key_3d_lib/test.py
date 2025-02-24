# Necessary Libraries
import vtk
import vtkmodules.vtkInteractionStyle
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkFiltersSources import vtkSphereSource
from vtkmodules.vtkFiltersSources import vtkConeSource
from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderWindow, vtkRenderer

# Create an instance of named colors in order to select colors for bg and object
colors = vtkNamedColors()

# Create and instance of VTK cone (pre-coded module)
def newCone(height, radius, res):
    cone = vtkConeSource()
    cone.SetHeight(height)
    cone.SetRadius(radius)
    cone.SetResolution(res)
    return cone

cone = newCone(2.0, 1.0, 10) # New cone instance

# Mapper process object needed to map polygonal data into graphics
# Output of the cone source connects to the input of the mapper
# 3D Graphics - Output shape should connect to the input of the mapper to be displayed

coneMapper = vtkPolyDataMapper()
coneMapper.SetInputConnection(cone.GetOutputPort())

# Actor renders the mapper's graphic primitives
# Also includes an internal transformation matrix

coneActor = vtkActor()
coneActor.SetMapper(coneMapper)
coneActor.GetProperty().SetColor(colors.GetColor3d('MistyRose'))

# Create the renderer and assign the actors to it
# Renderer - A part of the window on the screen and draws actors it has (also set bg here)

ren1 = vtkRenderer()
ren1.AddActor(coneActor)
ren1.SetBackground(colors.GetColor3d('MidnightBlue'))

# Create the renderer in a new window
renWin = vtkRenderWindow()
renWin.AddRenderer(ren1)
renWin.SetSize(300,300)

# interactor - should hold the window
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renWin)

renWin.Render()
interactor.Start()