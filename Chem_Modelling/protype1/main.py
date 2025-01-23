from input import Input

start = input("Y/N: ")

while start == "Y":
    info = input("Enter the molecule: ")
    hello = Input(info)
    hello.readInput
    hello.returnInput
    start = input("Y/N: ")
    if start == "N":
        break