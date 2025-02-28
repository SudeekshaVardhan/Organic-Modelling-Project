import datascraper
import pubchempy as pcp

inp = input("Start? (Y/N):")
inp = inp.lower()

newData = datascraper.Datascraper

index = 0

print(pcp.NotFoundError)

while inp != "n":

    newData.takeInput(newData)
    newData.datascraping(newData)
    print(newData.getCoords(newData))

    newData.molecule = None
    inp = input("Continue? (Y/N)")

    index += 1



    
