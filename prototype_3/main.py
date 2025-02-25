import datascraper


inp = input("Start? (Y/N):")
inp = inp.lower()

newData = datascraper.Datascraper

index = 0

while inp != "n":

    newData.takeInput(newData)
    newData.datascraping(newData)
    print(newData.getCoords(newData))


    newData.molecule = None
    inp = input("Continue? (Y/N)")

    index += 1

# newData.returnMF("ammonium thioglycolate")
# print(newData.getCoords(newData))


    
