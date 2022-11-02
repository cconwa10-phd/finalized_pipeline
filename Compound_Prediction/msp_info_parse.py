import pandas as pd
import requests
import time


def get_smiles_from_inchikey(inchikey):
    time.sleep(0.5)
    r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/'
                     f'compound/inchikey/{inchikey}/property/'
                     f'CanonicalSMILES,InChI/JSON')
    code = r.status_code
    print(code)
    if code != 200:
        csmiles = "Nah"
    else:
        r = r.json()
        csmiles = r['PropertyTable']['Properties'][0]['CanonicalSMILES']
    return csmiles


def gather_target(mspFile):
    nam = []
    inchk = []
    smile = []
    for line in open(mspFile):
        line = line.strip()
        if len(line) > 0:
            if ":" in line:
                key, value = line.split(":", 1)
                if key.strip() == "Name":
                    print(value.strip())
                    nam.append(value.strip())
                elif key.strip() == "InChIKey":
                    inchk.append(value.strip())
                    smile.append(get_smiles_from_inchikey(value.strip()))
                else:
                    pass
            else:
                pass
        else:
            pass
    df = pd.DataFrame(list(zip(nam, inchk, smile)), columns=["name", "inchikey", "smile"])
    return df


hold = 0
gather_target("/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/cmpd_test_software.msp").to_csv("results.csv")
