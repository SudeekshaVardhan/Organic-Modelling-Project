# UI Libraries
import tkinter as tk
from tkinter import PhotoImage
from PIL import Image, ImageTk
# Import the classes that get and process information
import datascraper # Imports MolSearch and MolModelling
import networkSearch # Imports NetSearch and NetModelling
# To access files saved in directory
import os

class Front:
    '''
    The main purpose of this class is to set up the UI. It uses Tkinter, an existing library.
    The Front class was created to satisfy Success Criteria #1: creating a easily accessible interface
    for teachers and students to use, especially if they are not familiar with terminals.
    
    '''
    def __init__(self):
        # Create a window object (customized - set color and size)
        self.root = tk.Tk()
        self.root.geometry('800x500') # Size of UI Window (not full screen)
        self.root.title('Prolycule') # Title label
        self.root.attributes('-topmost',True) # Ensure it appears on top of other windows

        # Nav frames - allows the user to flip between several frames
        self.main_frame = tk.Frame(self.root)
        self.new_wind = tk.Frame(self.root)
        self.main_frame.pack(fill="both", expand=True)

        # Other important variables
        self.input = tk.StringVar()
        self.data = datascraper.MolModelling()
        self.nets = networkSearch.NetModel()
    
        self.customWind()

        self.root.mainloop()

    def customWind(self):
        '''
        START WINDOW
        '''
        # Adds specific attributs to first frame - text, font, other stylistic elements
        label = tk.Label(self.main_frame, text='PROLYCULE', font=('Courier', 64))
        label.pack(padx= 20, pady=20)

        # Button to start process
        button1 = tk.Button(self.main_frame, text='START', font=('Courier', 16), command=self.newWind)
        button1.pack(padx=10, pady=10)

        # Exit button
        button2 = tk.Button(self.main_frame, text='CLOSE', font=('Courier', 16), command=self.root.destroy)
        button2.pack(padx=10,pady=10)

    def newWind(self):
        '''
        Select from network or molecule
        - Network EX: Buckminsterfullerene, graphene
        - Molecule EX: CO2, ethanol, acetate (ion)
        '''
        # Clear starting dialogue to page 1 (selections)
        for widget in self.main_frame.winfo_children():
            widget.destroy() # Sets each page up as a 'widget' so that when you move through them, you don't destroy

        # Prompt dialogue
        q1 = tk.Label(self.main_frame, text="Are you trying to access a network solid or a molecule?", font=('Courier', 16))
        q1.pack(padx=20,pady=20)

        # Buttons that allow user to choose between network and molecule
        button3 = tk.Button(self.main_frame, text='Network', font=('Courier', 16), command=self.networkChecker)
        button4 = tk.Button(self.main_frame, text='Molecule', font=('Courier', 16), command=self.molChecker)
        button3.pack(padx=10,pady=10)
        button4.pack(padx=10,pady=10)

        # Button to exit program fully
        button5 = tk.Button(self.main_frame, text="Exit", font=("Courier", 12), command=self.root.destroy) # Kills entire program, not just root
        button5.pack(padx=10, pady=100)

    def molChecker(self):
        '''Input for molecular compounds
        Takes the name/formula, outputs a 2D image and a 3D strucure of the molecule

        
        '''
        for widget in self.main_frame.winfo_children():
            widget.destroy() # Widget to enter the molecule (if user chooses the molecule on newWind, navigate here)
        
        # Prompt dialogue
        quer1 = tk.Label(self.main_frame, text= "Enter name or formula of molecule: ", font=('Courier', 16))
        quer1.pack()

        # Enter the mol formula or name (both are applicable in PubChempy)
        myentry = tk.Entry(self.main_frame, textvariable= self.input, font=('Courier', 16), width=20)
        myentry.pack()

        # Submit button should access the datascraping method of the Datascraper class and return image
        def submit():
            self.input.set(myentry.get())
            print(self.input.get())
            self.molViewer()

        # Submission button
        button6 = tk.Button(self.main_frame, text="Submit", font=('Courier', 16), command=submit)
        button6.pack(padx=10,pady=10)

        # Return to selection page (if user wants to navigate back to check a network)
        button7 = tk.Button(self.main_frame, text="Back", font=("Courier", 12),command=self.newWind)
        button7.pack(padx=100,pady=150)

    def networkChecker(self):
        '''
        Input for network solids
        This takes the name or formula of a network solid, passes it through the networkSearch class and PubChempy
        to get the 2D structure and the 3D strucutre.
        
        Builds off the networkSearch/netModelling class
        '''

        for widget in self.main_frame.winfo_children():
            widget.destroy() # Navigate here if user selects the network button on the newWind widget

        # Prompt dialogue
        quer2 = tk.Label(self.main_frame, text= "Enter name or formula of solid or allotrope: ", font=('Courier', 16))
        quer2.pack()

        # Enter formula or name (both are applicable in the COD)
        myentry = tk.Entry(self.main_frame, textvariable=self.input, font=('Courier', 16), width=20)
        myentry.pack()

        # Submit button should access the network database to get the dat from COD
        def submit():
            self.input.set(myentry.get())
            print(self.input.get())
            self.netViewer()

        button6 = tk.Button(self.main_frame, text="Submit", font=('Courier', 16), command=submit)
        button6.pack(padx=10,pady=10)

        # Return to selection page
        button7 = tk.Button(self.main_frame, text="Back", font=("Courier", 12),command=self.newWind)
        button7.pack(padx=100,pady=150)
        
    def molViewer(self):
        '''
        New frame to view molecule
        Displays the molecule image and the other buttons for additional information
        '''
        for widget in self.main_frame.winfo_children():
            widget.destroy() # Navigate here agter the submit button is pressed on molChecker

        # Creates instance of Datascraper class to access the image-generating datascraper class
        self.data.datascraping(self.input.get())

        # Get the correct file path for the mol image
        script_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(script_dir, "mol.png")

        if not os.path.exists(file_path): # If image path is invalid (debugging node)
            print("Error: Image file not found")
            return
        # Display 2D rendering created by the Datascraper class
        try:
            # Convert image file to tkinter readable format
            img = Image.open(file_path)
            img = img.resize((200, 200), Image.Resampling.LANCZOS)
            image = ImageTk.PhotoImage(img)
            
            # Display image as a label and add the label to tkinter
            im_label = tk.Label(self.main_frame, image=image)
            im_label.image = image
            self.image_ref = image
            im_label.pack(pady=20)

        # Error handling
        except FileNotFoundError:
            print("Error: mol.png not found")
        
        except Exception as e:
            print(f"Unexpected error occured: {e}")

        # User can select additional information about the molecule
        buttonA = tk.Button(self.main_frame, text="Other Names", font=('Courier', 12), command=self.returnMolAltNames)
        buttonB = tk.Button(self.main_frame, text="Properties", font=('Courier', 12), command=self.returnMolProperties)
        buttonA.pack(padx=10, pady=10)
        buttonB.pack(padx=10, pady=10)

        '''
        Key for molecules:
            Red = Oxygen
            Gray = Carbon
            White = Hydrogen
            Blue = Nitrogen
            Yellow = Sulfur
        '''
        label = tk.Label(self.main_frame, text="KEY", font=("Courier", 16))
        carb = tk.Label(self.main_frame, text="Carbon = Gray", font=("Courier", 12))
        oxy = tk.Label(self.main_frame, text="Oxygen = Red", font=("Courier", 12))
        hyd = tk.Label(self.main_frame, text="Hydrogen = White", font=("Courier", 12))
        nitro = tk.Label(self.main_frame, text="Nitro = Blue", font=("Courier", 12))
        sulfur = tk.Label(self.main_frame, text="Sulfur = Yellow", font=("Courier", 12))

        label.pack(padx=10,pady=10)
        carb.pack(padx=10,pady=10)
        oxy.pack(padx=10,pady=10)
        hyd.pack(padx=10,pady=10)
        nitro.pack(padx=10,pady=10)
        sulfur.pack(padx=10,pady=10)

        # Button to navigate back to molChecker if user wishes to input another molecule
        returnHome = tk.Button(self.main_frame, text="Home?", font=('Courier', 12), 
                               command=self.molChecker)
        returnHome.pack(padx=10,pady=10)

        # 3D Modelling functionality
        # Must be at the end so it doesn't interrupt other processes.
        self.data.renWin()      
          
    def netViewer(self):
        for widget in self.main_frame.winfo_children(): 
            widget.destroy() # Navigate here once submit is pressed in netChecker()

        # Title + additional information
        txt = tk.Label(self.main_frame, 
                      text= "Note: The image displays a single repeating unit of a network.", 
                      font=("Courier", 10))
        txt.pack(padx=10, pady=10)

        txt2 = tk.Label(self.main_frame, 
                      text= "Refer to the structure to observe the network as a whole. Individual atoms are not shown", 
                      font=("Courier", 10))
        txt2.pack(padx=10, pady=10)


        # Instantiates the command to get the single repeating unit of the network
        self.nets.getCIF(self.input.get())

        # Get the correct file path for the mol image
        script_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(script_dir, "network.png")

        # Display 2D rendering created by the Datascraper class
        try:
            # Convert image file to tkinter readable format
            img = Image.open(file_path)
            img = img.resize((200, 200), Image.Resampling.LANCZOS)
            image = ImageTk.PhotoImage(img)
            
            # Display image as a label
            im_label = tk.Label(self.main_frame, image=image)
            im_label.image = image
            im_label.pack(pady=20)

            # Refresh the image after each call of the data class
            self.root.update_idletasks()  
            self.root.update()

        # Error handling
        except FileNotFoundError: # If the image is not found (to indicate wrong path)
            print("Error: mol.png not found")
        
        except Exception as e: # Any additional error occured
            print(f"Unexpected error occured: {e}")
        
        # User can choose whether they want to look at other 
        buttonA = tk.Button(self.main_frame, text="Other Names", font=('Courier', 12), command=self.returnNetAltNames)
        buttonB = tk.Button(self.main_frame, text="Properties", font=('Courier', 12), command=self.returnNetProperties)
        buttonA.pack(padx=10, pady=10)
        buttonB.pack(padx=10, pady=10)


        # Navigates back to the home screen
        returnHome = tk.Button(self.main_frame, text="Home?", font=('Courier', 12), 
                               command=self.molChecker)
        returnHome.pack(padx=10,pady=10)

        # Call the method necessary to 3D model the network
        self.nets.netRenWin()

    def returnMolAltNames(self):
        for widget in self.main_frame.winfo_children():
            widget.destroy() # Navigate here once Other Names is pressed on molViewer window
        
        # Title
        text = f"Alternate Names for {self.input.get()}"
        molTitle = tk.Label(self.main_frame, text=text, font=("Courier", 16))
        molTitle.pack(padx=10,pady=10)

        # Create a listbox for users to scroll through all the names 
        # This was chosen due to its adaptable size (ideal for varying datasets of alternate names for different compounds)
        box1 = tk.Listbox(self.main_frame, width=70) 
        box1.pack()

        # Add names to the listbox
        list = self.data.returnNames(self.input.get())
        for data in list:
            box1.insert(tk.END, data)

        # Return to molViewer
        homeButton = tk.Button(self.main_frame, text="Back", font=('Courier', 12), command=self.molViewer)
        homeButton.pack(padx=10, pady=10)

    def returnMolProperties(self):
        '''
        Return certain properties for the entered compound. Refer to the getProperties() method in the data class
        to get a full list of everything shown.

        Note: if the network or compound does not contain data, an error will be returned. In addition, not all values
        may be filled for every compound.    
        '''
        for widget in self.main_frame.winfo_children():
            widget.destroy() # Navigate here from molViewer once the Properties button is selected

        # Title
        text = f"Properties for {self.input.get()}"
        molTitle = tk.Label(self.main_frame, text=text, font=("Courier", 16))
        molTitle.pack(padx=10, pady=10)

        info = str(self.input.get()).lower()
        properties = [
            f"Log P Value: {self.data.getProperties(info, 'logP')}",
            f"Compound Charge: {self.data.getProperties(info, 'charge')}",
            f"IUPAC Name: {self.data.getProperties(info, 'iupac')}",
            f"Boiling Point: {self.data.getProperties(info, 'boiling_point')}",
            f"Melting Point: {self.data.getProperties(info, 'melting_point')}",
            f"Specific Heat Capacity: {self.data.getProperties(info, 'specific_heat_capacity')}",
            f"Enthalpy at 298 K (): {self.data.getProperties(info, 'enthalpy')}",
            f"Entropy at 298 K (): {self.data.getProperties(info, 'entropy')}",
            f"Gibbs Free Energy at 298 K (): {self.data.getProperties(info, 'gibbs_free_energy')}",  
        ]

        for val in properties: # Loop through the array and add value to UI screen
            label = tk.Label(self.main_frame, text=val, font=("Courier", 12))
            label.pack()

        # Return to molViewer
        homeButton = tk.Button(self.main_frame, text="Back", font = ("Courier", 12), command=self.molViewer)
        homeButton.pack(padx=10, pady=10)
    
    def returnNetAltNames(self):
        '''
        Return any other names that may exist for a network solid. These range from common names to complex chemical names.
        Like with molecules, the names are sorted from shortest to longest.
        '''
        for widget in self.main_frame.winfo_children():
            widget.destroy() # Navigate here from netViewer if Names button is selected
        
        # Title
        text = f"Alternate Names for {self.input.get()}"
        molTitle = tk.Label(self.main_frame, text=text, font=("Courier", 16))
        molTitle.pack(padx=10,pady=10)

        # Create a listbox for users to scroll through all the names 
        # This was chosen due to its adaptable size (ideal for varying datasets of alternate names for different compounds)
        box1 = tk.Listbox(self.main_frame, width=70) 
        box1.pack()
        
        # Add names to the listbox
        list = self.data.returnNames(self.input.get())
        for data in list:
            box1.insert(tk.END, data)
        
        # Navigate back to the netViewer page (if user wishes to go back further they can from there)
        homeButton = tk.Button(self.main_frame, text="Back", font=("Courier", 12), command=self.netViewer)
        homeButton.pack(padx=10, pady=10)

    def returnNetProperties(self):
        '''
        Return certain properties for the entered compound. Refer to the getProperties() method in the data class
        to get a full list of everything shown.

        Note: if the network or compound does not contain data, an error will be returned. In addition, not all values
        may be filled for every compound.    
        '''
        for widget in self.main_frame.winfo_children():
            widget.destroy() # Navigate here if the user wishes to view properties for their network
        text = f"Properties for {self.input.get()}" # Widget title
        molTitle = tk.Label(self.main_frame, text=text, font=("Courier", 16))
        molTitle.pack(padx=10, pady=10)

        info = str(self.input.get()).lower() # Ensure input is stripped for easier processing

        # List of properties displayed in the call
        properties = [
            f"Log P Value: {self.data.getProperties(info, 'logP')}",
            f"Compound Charge: {self.data.getProperties(info, 'charge')}",
            f"IUPAC Name: {self.data.getProperties(info, 'iupac')}",
            f"Boiling Point: {self.data.getProperties(info, 'boiling_point')}",
            f"Melting Point: {self.data.getProperties(info, 'melting_point')}",
            f"Specific Heat Capacity: {self.data.getProperties(info, 'specific_heat_capacity')}",
            f"Enthalpy at 298 K (): {self.data.getProperties(info, 'enthalpy')}",
            f"Entropy at 298 K (): {self.data.getProperties(info, 'entropy')}",
            f"Gibbs Free Energy at 298 K (): {self.data.getProperties(info, 'gibbs_free_energy')}",  
        ]

        for val in properties: # Loop through the array and add value to UI screen
            label = tk.Label(self.main_frame, text=val, font=("Courier", 12))
            label.pack()

        # Return to the netViewer screen
        homeButton = tk.Button(self.main_frame, text="Back", font = ("Courier", 12), command=self.netViewer)
        homeButton.pack(padx=10, pady=10)

front1 = Front()