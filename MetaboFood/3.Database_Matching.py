import datetime
import requests
import pandas as pd
import time

#Need database csv file, need hmdb file, need target spectra file to search against

def get_smiles_from_name(name): #matches on name
    r = requests.get(
        f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES,InChI,IUPACName,InChIKey/JSON')
    code = r.status_code
    print(code)
    if code != 200:
        csmiles = "Nah"
        cinchi = "Nah"
        iupac = "Nah"
        inchikey = "Nah"
    else:
        r = r.json()
        csmiles = r['PropertyTable']['Properties'][0]['CanonicalSMILES']
        cinchi = r['PropertyTable']['Properties'][0]['InChI']
        iupac = r['PropertyTable']['Properties'][0]['IUPACName']
        inchikey = r['PropertyTable']['Properties'][0]['InChIKey']
        #title = r['PropertyTable']['Properties'][0]['Title']
    return csmiles, cinchi, iupac, inchikey


def inchi_list_USDA(file):
    inchi_l = []
    smile_l = []
    iupac_l = []
    inchiK_l = []
    data = pd.read_csv(file, encoding = "ISO-8859-1")
    names = data["compound"].values
    for i, name in enumerate(names):
        s, ih, iu, ik = get_smiles_from_name(name)
        inchi_l.append(ih)
        smile_l.append(s)
        iupac_l.append(iu)
        inchiK_l.append(ik)
        time.sleep(0.5)
        print(datetime.datetime.now().strftime("%H:%M:%S"))
    data["name"] = names
    data["InChI"] = inchi_l
    data["SMILE"] = smile_l
    data["IUPAC"] = iupac_l
    data["InChIKey_use"] = inchiK_l
    return data

def get_info_from_inchikey(inchikey):
    r = requests.get(
        f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/CanonicalSMILES,InChI,IUPACName/JSON')
    code = r.status_code
    print(code)
    if code != 200:
        csmiles = "Nah"
        cinchi = "Nah"
        iupac = "Nah"
    else:
        r = r.json()
        if not (r['PropertyTable']['Properties'][0]['CanonicalSMILES'] is None):
            csmiles = r['PropertyTable']['Properties'][0]['CanonicalSMILES']
        else:
            csmiles = "Nah"
        if not(r['PropertyTable']['Properties'][0]['InChI'] is None):
            cinchi = r['PropertyTable']['Properties'][0]['InChI']
        else:
            cinchi = "Nah"
        # if not (r['PropertyTable']['Properties'][0]['IUPACName'] is None):
        #     iupac = r['PropertyTable']['Properties'][0]['IUPACName']
        # else:
        #     iupac = "Nah"
        #title = r['PropertyTable']['Properties'][0]['Title']
    return csmiles, cinchi

def inchi_list_new(file):
    inchi_l = []
    smile_l = []
    data = pd.read_csv(file)
    inchikeys = data["InChIKey"].values
    for i, inchikey in enumerate(inchikeys):
        s, i = get_info_from_inchikey(inchikey)
        inchi_l.append(i)
        smile_l.append(s)
        time.sleep(0.5)
        print(datetime.datetime.now().strftime("%H:%M:%S"))
    data["InChiKey2"] = inchikeys
    data["InChI2"] = inchi_l
    data["SMILES2"] = smile_l
    return data

def hmdb_xml_match_merge_name(target_df, ref_df, filename):
    target_left = pd.read_csv(target_df, dtype=str, encoding = "ISO-8859-1")
    ref_right = pd.read_csv(ref_df, dtype=str, encoding = "ISO-8859-1")
    print(target_left.head(), ref_right.head())
    hmdb_matches = pd.merge(target_left, ref_right, how="left", on=["compound"])
    hmdb_matches.to_csv(filename)
    print(hmdb_matches.head())
    return hmdb_matches

def hmdb_xml_match_merge_inchikey(target_df, ref_df, filename):
    target_left = pd.read_csv(target_df, dtype=str, encoding = "ISO-8859-1")
    ref_right = pd.read_csv(ref_df, dtype=str, encoding = "ISO-8859-1")
    print(target_left.head(), ref_right.head())
    hmdb_matches = pd.merge(target_left, ref_right, how="left", on=["inchikey"])
    hmdb_matches.to_csv(filename)
    print(hmdb_matches.head())
    return hmdb_matches

def inch_loop(data_inchks, msp_inchk):
    found = False
    for i, inchik in enumerate(data_inchks):
        if inchik.strip() == msp_inchk:
            found = True
    return found

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def assign_spectra(filename, msp):
    mspFile = ""
    mzval = []
    intval = []
    meta_data = {}
    data = pd.read_csv(filename, encoding = "ISO-8859-1")
    inchiks = data["InChIKey_use"].values

    for line in open(msp):
        line = line.strip()
        if len(line) > 0:
            if ":" in line:
                key, value = line.split(":", 1)
                meta_data[key] = value
            else:
                mz, intensity = line.split(' ', 1)
                mzval.append(float(mz.strip()))
                if isfloat(intensity):
                    intval.append(float(intensity.strip()))
                else:
                    inten, extra = intensity.split(" ", 1)
                    intval.append(float(inten.strip()))
        else:
            if meta_data.get("InChIKey") is not None:
                found = inch_loop(inchiks, meta_data.get("InChIKey").strip())
                if found:
                    for kMeta, vMeta in meta_data.items():
                        mspFile += str(kMeta) + ": " + str(vMeta) + "\n"
                    for kMZ, vInt in zip(mzval, intval):
                        mspFile += str(kMZ) + "\t" + str(vInt) + "\n"
                    mspFile += "\n"
                    print(meta_data)
                else:
                    pass
            else:
                pass
            mzval = []
            intval = []
            meta_data = {}
    print(mspFile)
    return mspFile

def write_out_msp_USDA(mspFile, file):
    out_file = open(file, "w")
    out_file.write(mspFile)
    out_file.close()
    print("Outfile: ", file, " has been saved")


def main():
    pass


if __name__ == '__main__':
    main()
