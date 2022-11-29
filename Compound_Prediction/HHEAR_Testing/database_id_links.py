import pandas as pd
import numpy as np
import requests
import time
from collections import defaultdict

def get_info_from_id(id_):
    r = requests.get(
        f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{id_}/property/CanonicalSMILES,InChI,IUPACName,InChIKey/JSON')
    code = r.status_code
    print(code)
    if code != 200:
        csmiles = "Nah"
        cinchi = "Nah"
    else:
        r = r.json()
        if not (r['PropertyTable']['Properties'][0]['CanonicalSMILES'] is None):
            csmiles = r['PropertyTable']['Properties'][0]['CanonicalSMILES']
        else:
            csmiles = "Nah"
        if not(r['PropertyTable']['Properties'][0]['InChIKey'] is None):
            cinchi = r['PropertyTable']['Properties'][0]['InChIKey']
        else:
            cinchi = "Nah"
    return csmiles, cinchi

def hmdb_xml_match_merge(target_left, ref_df, filname):
    #target_left = pd.read_csv(target_df, dtype=str)
    ref_right = pd.read_csv(ref_df, dtype=str)
    print(target_left.head(), ref_right.head())
    hmdb_matches = pd.merge(target_left, ref_right, how="left", on=["HMDB"])
    hmdb_matches.to_csv(filname)
    print(hmdb_matches.head())
    return hmdb_matches

def inchi_key_smiles(mspFile, ref_data):
    #8===D -- -- -- david did this
    metaData = defaultdict(list)
    for line in open(mspFile):
        line = line.strip()
        if len(line) > 0 and ":" in line:
            key, value = line.split(":", 1)
            if key.strip() == "Name":
                metaData[key].append(value.strip())
            elif key.strip() == "DATABASE_ID" or key.strip() == "PUBCHEM_COMPOUND_CID":
                time.sleep(0.5)
                print(value.strip())
                s, i = get_info_from_id(value.strip())
                metaData["HMDB"].append(value.strip())
                metaData["SMILES"].append(s)
                metaData["InChIKey"].append(i)
        else:
            pass
    df = pd.DataFrame(metaData)
    inchi_col = df.pop("HMDB")
    df.insert(0,"HMDB", inchi_col)
    df1 = hmdb_xml_match_merge(df, ref_data, "target_compounds2.csv")
    return df1



def main():
    inchi_key_smiles("/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/hhear2.msp", "/Users/ciaraconway/Documents/all_databases/Compounds/full_hmdb_xml_parsed_inchikey.csv")


if __name__ == '__main__':
    main()
