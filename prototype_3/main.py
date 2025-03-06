import datascraper
import pubchempy as pcp

inp = input("Start? (Y/N):")
inp = inp.lower()

newData = datascraper.Datascraper
newRen = datascraper.Modelling(newData)


print(pcp.NotFoundError)

while inp != "n":

    newData.takeInput(newData)
    newData.datascraping(newData)

    newRen.renWin()

    newData.molecule = None
    inp = input("Continue? (Y/N)")




    
