import requests
import os

ACCESS_TOKEN = os.getenv("ZENODO_ACCESS_TOKEN")

headers = {"Content-Type": "application/json"}
params = {'access_token': ACCESS_TOKEN}

r = requests.get(
    'https://zenodo.org/api/deposit/depositions',
    params=params,
    #json={},
    #headers=headers,
)

print(r.status_code)
print(r.json())

