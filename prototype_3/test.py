import requests
import json

class NetworkSearch:
    def __init__(self):
        self.network = None
        self.CODID = None
        
    def getCODID(self, name):
        self.network = name
        url = f"https://www.crystallography.net/cod/result?text={self.network}&format=json"
        response = requests.get(url)

        if response.status_code == 200:  # If the API request succeeds
            try:
                data = response.json()

                if isinstance(data, list) and data:  # Ensure it's a list with results
                    self.CODID = data[0].get("file", "No COD ID found")
                    return self.CODID
                else:
                    print("No valid COD entries found.")
                    return None

            except json.JSONDecodeError:
                print("Error: Failed to decode JSON response.")
                return None
        else:
            print(f"Error fetching data. Status code: {response.status_code}")
            return None

# Example Usage
search = NetworkSearch()
compound_name = input("Enter a compound name: ")  # User enters compound name
cod_id = search.getCODID(compound_name)

if cod_id:
    print(f"First COD ID found for {compound_name}: {cod_id}")
else:
    print("No COD ID found.")