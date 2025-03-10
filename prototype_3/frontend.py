import tkinter as tk
import datascraper

class Front:
    
    root = tk.Tk()
    def __init__(self):
        pass
        
    def createWind (self):
        # Create a window object
        self.root.geometry('800x500')
        self.root.title('Prolycule')

    def inputStuff(self):
        label = tk.Label(self.root, text='PROLYCULE', font=('Courier', 64))
        label.pack(padx= 20, pady=20)

        # Create box for text + entry (not used)
        '''
        text = tk.Text(root, font=('Courier', 16))
        text.pack()
        myentry = tk.Entry(root)
        myentry.pack()
        '''

        # Button
        def buttonPressed():
            global buttonCall
            buttonCall = not buttonCall

        buttonCall = False

        button = tk.Button(root, text='START', font=('Courier', 16), command=buttonPressed())
        button.pack(padx=10, pady=10)

        if buttonCall:
            print("I'm being pressed!")

            
        # Last 
        self.root.mainloop()