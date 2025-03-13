import tkinter as tk
from tkinter import PhotoImage
from PIL import Image, ImageTk
import datascraper
import os

class Front:
    
    def __init__(self):
        # Create a window object (customized)
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
    
        self.customWind()

        self.root.mainloop()


    def customWind(self):
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
        - Molecule EX: 
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
        button7 = tk.Button(self.main_frame, text="Back", font=("Courier", 16),command=self.newWind)
        button7.pack(pady=100)

    def networkChecker(self):
        '''Input for network solids'''
        print("This is network checker")

    def molViewer(self):
        '''New frame to view molecule. Add additional data'''
        for widget in self.main_frame.winfo_children():
            widget.destroy()

        # Creates instance of Datascraper class to access the image-generating datascraper class
        gen = datascraper.Datascraper()
        gen.datascraping(self.input.get())

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
            im_label.image = image
            im_label.pack(pady=20)

        except FileNotFoundError:
            print("Error: mol.png not found")
        
        except Exception as e:
            print(f"Unexpected error occured: {e}")
        
        buttonA = tk.Button(self.main_frame, text="")
        buttonB = tk.Button(self.main_frame, text="")
        buttonC = tk.Button(self.main_frame, text="")
        buttonA.pack()
        buttonB.pack()
        buttonC.pack()