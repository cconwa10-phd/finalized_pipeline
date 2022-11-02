import os
import subprocess
import datetime
import requests
import pandas as pd
import time


def get_smiles_from_name(name):
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
    names = data["NutrDesc"].values
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


def hmdb_xml_match_merge_name(target_df, ref_df, filename):
    target_left = pd.read_csv(target_df, dtype=str, encoding = "ISO-8859-1")
    ref_right = pd.read_csv(ref_df, dtype=str, encoding = "ISO-8859-1")
    print(target_left.head(), ref_right.head())
    hmdb_matches = pd.merge(target_left, ref_right, how="left", on=["InChIKey_use"])
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
    inchiks = data["InChIKey"].values

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
    # file = "/Users/ciaraconway/Documents/Kay_Lab/Taxonomy_FullDataSet_2021-10-08_Inchi_smiles.csv"
    # ref_file = "/Users/ciaraconway/Documents/HHEAR_BIO/full_hmdb_xml_parsed_name.csv"
    # nist_file = "/Users/ciaraconway/Documents/Kay_Lab/Taxonomy_FullDataSet_2021-10-08_Inchi_smiles_nist.csv"
    # hmdb_file = "/Users/ciaraconway/Documents/Kay_Lab/Taxonomy_FullDataSet_2021-10-08_Inchi_smiles_hmdb.csv"
    # msp_file = "/Users/ciaraconway/Documents/Kay_Matching/hr_msms_nist.MSP"
    # new_msp_file = "/Users/ciaraconway/Documents/Kay_Lab/Taxonomy_FullDataSet_2021-10-08_Inchi_smiles.msp"

    #file = "/Users/ciaraconway/Documents/Kay_Matching/NUTR_DEF_pc.csv"
    file = "/Users/ciaraconway/Documents/Kay_Matching/nutrient_new.csv"
    # ref_file = "/Users/ciaraconway/Documents/HHEAR_BIO/full_hmdb_xml_parsed_name.csv"
    # nist_file = "/Users/ciaraconway/Documents/Kay_Matching/nutrient_nist_new.csv"
    # hmdb_file = "/Users/ciaraconway/Documents/Kay_Matching/nutrient_hmdb_new.csv"
    hhear_file = "/Users/ciaraconway/Desktop/inchikey_search_HHEAR.csv"
    msp_file = "/Users/ciaraconway/Documents/Kay_Matching/hr_msms_nist.MSP"
    msp_file_mona_lcmsmspn = "/Users/ciaraconway/Desktop/lcmsmspn.msp"
    new_msp_file = "/Users/ciaraconway/Desktop/LCMSMS_Exposome_MoNA.msp"

    #inchi_list_USDA(file).to_csv(nist_file)
    #hmdb_xml_match_merge_name(nist_file, ref_file, hmdb_file)
    new_msp = assign_spectra(hhear_file, msp_file_mona_lcmsmspn)
    write_out_msp_USDA(new_msp, new_msp_file)








if __name__ == '__main__':
    main()

