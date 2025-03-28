# the 
import tkinter as tk
from tkinter import PhotoImage
from PIL import Image, ImageTk
# Import the classes that get and process information
import datascraper # Imports MolSearch and MolModelling
import networkSearch # Imports
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
        self.root.geometry('800x500')
        self.root.title('Prolycule')
        self.root.attributes('-topmost',True)

        # Nav frames
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
        # Adds specific attributs to first frame
        label = tk.Label(self.main_frame, text='PROLYCULE', font=('Courier', 64))
        label.pack(padx= 20, pady=20)

        # Button to start process
        button1 = tk.Button(self.main_frame, text='START', font=('Courier', 16), command=self.newWind)
        button1.pack(padx=10, pady=10)

        button2 = tk.Button(self.main_frame, text='CLOSE', font=('Courier', 16), command=self.root.destroy)
        button2.pack(padx=10,pady=10)

        # Images

    def newWind(self):
        '''
        Select from network or molecule
        - Network EX: Buckminsterfullerene, graphene
        - Molecule EX: CO2, ethanol, acetate (ion)
        '''
        # Clear starting dialogue to page 1 (selections)
        for widget in self.main_frame.winfo_children():
            widget.destroy()

        q1 = tk.Label(self.main_frame, text="Are you trying to access a network solid or a molecule?", font=('Courier', 16))
        q1.pack(padx=20,pady=20)

        button3 = tk.Button(self.main_frame, text='Network', font=('Courier', 16), command=self.networkChecker)
        button4 = tk.Button(self.main_frame, text='Molecule', font=('Courier', 16), command=self.molChecker)

        button3.pack(padx=10,pady=10)
        button4.pack(padx=10,pady=10)

        # Button to exit program fully
        button5 = tk.Button(self.main_frame, text="Exit", font=("Courier", 12), command=self.root.destroy)
        button5.pack(padx=10, pady=100)

    def molChecker(self):
        '''Input for molecular compounds'''
        for widget in self.main_frame.winfo_children():
            widget.destroy()
        
        # Mol query: enter name or MF
        quer1 = tk.Label(self.main_frame, text= "Enter name or formula of molecule: ", font=('Courier', 16))
        quer1.pack()

        myentry = tk.Entry(self.main_frame, textvariable= self.input, font=('Courier', 16), width=20)
        myentry.pack()

        # Submit button should access the datascraping method of the Datascraper class and return image
        def submit():
            self.input.set(myentry.get())
            print(self.input.get())
            self.molViewer()

        button6 = tk.Button(self.main_frame, text="Submit", font=('Courier', 16), command=submit)
        button6.pack(padx=10,pady=10)

        # Return to selection page
        button7 = tk.Button(self.main_frame, text="Back", font=("Courier", 12),command=self.newWind)
        button7.pack(padx=100,pady=150)

    def networkChecker(self):
        '''
        Input for network solids
        This takes the name or formula of a network solid, passes it through the networkSearch class and PubChempy
        to get the 2D structure and the 3D strucutre.
        
        '''

        for widget in self.main_frame.winfo_children():
            widget.destroy()

        quer2 = tk.Label(self.main_frame, text= "Enter name or formula of solid or allotrope: ", font=('Courier', 16))
        quer2.pack()

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
        '''New frame to view molecule. Add additional data'''
        for widget in self.main_frame.winfo_children():
            widget.destroy()

        # Creates instance of Datascraper class to access the image-generating datascraper class
        self.data.datascraping(self.input.get())

        # Get the correct file path for the mol image
        script_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(script_dir, "mol.png")

        # Display 2D rendering created by the Datascraper class
        try:
            # Convert image file to tkinter readable format
            img = Image.open(file_path)
            img = img.resize((200, 200), Image.Resampling.LANCZOS)
            image = ImageTk.PhotoImage(img)
            
            # Display image as a label
            im_label = tk.Label(self.main_frame, image=image)
            im_label.pack(pady=20)

        # Debugging stuff
        except FileNotFoundError:
            print("Error: mol.png not found")
        
        except Exception as e:
            print(f"Unexpected error occured: {e}")

        buttonA = tk.Button(self.main_frame, text="Other Names", font=('Courier', 12), command=self.returnAltNames)
        buttonB = tk.Button(self.main_frame, text="Properties", font=('Courier', 12), command=self.returnProperties)

        buttonA.pack(padx=10, pady=10)
        buttonB.pack(padx=10, pady=10)

        returnHome = tk.Button(self.main_frame, text="Home?", font=('Courier', 12), 
                               command=self.molChecker)
        returnHome.pack(padx=10,pady=10)

        # 3D Modelling functionality
        # Must be at the end so it doesn't interrupt other processes.
        self.data.renWin()      
          
    def netViewer(self):
        for widget in self.main_frame.winfo_children(): # Switch to the next page of the program
            widget.destroy()

        # Instantiates the 
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
            im_label.pack(pady=20)

        # Debugging stuff
        except FileNotFoundError:
            print("Error: mol.png not found")
        
        except Exception as e:
            print(f"Unexpected error occured: {e}")
        
        buttonA = tk.Button(self.main_frame, text="Other Names", font=('Courier', 12), command=self.returnAltNames)
        buttonB = tk.Button(self.main_frame, text="Properties", font=('Courier', 12), command=self.returnProperties)

        buttonA.pack(padx=10, pady=10)
        buttonB.pack(padx=10, pady=10)


        # Navigates back to the home screen
        returnHome = tk.Button(self.main_frame, text="Home?", font=('Courier', 12), 
                               command=self.molChecker)
        returnHome.pack(padx=10,pady=10)

        # Call the method necessary to 3D model the class
        self.nets.netRenWin()

    def returnMolAltNames(self):
        for widget in self.main_frame.winfo_children():
            widget.destroy()
        
        text = f"Alternate Names for {self.input.get()}"
        molTitle = tk.Label(self.main_frame, text=text, font=("Courier", 16))
        molTitle.pack(padx=10,pady=10)

        box1 = tk.Listbox(self.main_frame, width=50) 
        box1.pack()
            
        list = self.data.returnNames(self.input.get())
        for data in list:
            box1.insert(tk.END, data)

        homeButton = tk.Button(self.main_frame, text="Back", command=self.molViewer)
        homeButton.pack(padx=10, pady=10)

    def returnMolProperties(self):
        for widget in self.main_frame.winfo_children():
            widget.destroy()
        text = f"Properties for {self.input.get()}"
        molTitle = tk.Label(self.main_frame, text=text, font=("Courier", 16))
        molTitle.pack(padx=10, pady=10)

        properties = []
        label = tk.Label(self.main_frame, )

        homeButton = tk.Button(self.main_frame, text="Back", command=self.newWind)
        homeButton.pack(padx=10, pady=10)

front1 = Front()