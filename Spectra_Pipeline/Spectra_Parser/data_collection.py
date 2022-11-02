import requests
headers = {
    'Accept': 'text/msp',
}
response = requests.get('https://mona.fiehnlab.ucdavis.edu/rest/spectra/search?query=tags.text=="LC-MS" and metaData=q=\'name=="ms level" and value=="MS2"\'', headers=headers)
