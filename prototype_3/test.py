import requests
from mp_api import MPRester

search_url = "https://materialsproject.org/rest/v2/query"
query = {
    "criteria" : {"pretty_formula": "NaCl"},  # Change this to a simple formula you know exists
    "properties" : ["task_id", "pretty_formula"]
}

header = {"X-API-KEY" : "your_new_api_key_here"}
response = requests.post(search_url, json=query, headers=header)

if response.status_code == 200:
    print("Query success:", response.json())
else:
    print("Error:", response.status_code, response.text)