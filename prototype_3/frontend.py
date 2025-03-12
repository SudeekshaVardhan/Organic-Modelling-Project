import tkinter as tk
import datascraper

class Front:
    
    def __init__(self):
        # Create a window object 
        self.root = tk.Tk()
        self.root.geometry('800x500')
        self.root.title('Prolycule')
    
        self.customWind()

        self.root.mainloop()


    def customWind(self):
        label = tk.Label(self.root, text='PROLYCULE', font=('Courier', 64))
        label.pack(padx= 20, pady=20)

        # Button to start process
        button1 = tk.Button(self.root, text='START', font=('Courier', 16), command=self.newWind)
        button1.pack(padx=10, pady=10)

        button2 = tk.Button(self.root, text='CLOSE', font=('Courier', 16), command=self.root.destroy)
        button2.pack(padx=10,pady=10)

        # Images

    def newWind(self):
        new_wind = tk.Toplevel()
        new_wind.geometry("800x500")
        text = tk.Text(new_wind, font=('Courier', 16))
        text.pack()
        myentry = tk.Entry(new_wind)
        myentry.pack()
        