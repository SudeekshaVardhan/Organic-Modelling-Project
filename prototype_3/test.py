import requests
from bs4 import BeautifulSoup

# Search for a compound on NIST
compound_name = "methane"
search_url = f"https://webbook.nist.gov/cgi/cbook.cgi?Name={compound_name}&Units=SI"

response = requests.get(search_url)
soup = BeautifulSoup(response, 'html.parser')

# Extract enthalpy data from the page
enthalpy_data = soup.find(text="Enthalpy of formation")
if enthalpy_data:
    print(enthalpy_data.find_next().text)
else:
    print("Enthalpy data not found.")