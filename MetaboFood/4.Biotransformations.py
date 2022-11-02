import os
import subprocess
import datetime
import requests
import pandas as pd
import time

###Need to have Biotransformer Jar File path specified to call###

#-s,--nsteps <Number of steps> The number of steps for the prediction. This option can be set by the user for the EC-based, CYP450, Phase II, and Environmental microbial biotransformers. The default value is 1.

#Can test connection to Biotransformer with this
# os.system("echo hello")
# os.system("cd biotransformer3.0jar/")
# os.system("ls")
# print(os.getcwd())
# subprocess.run(["cd", "biotransformer3.0jar/"])
# os.chdir("biotransformer3.0jar")
# subprocess.run(["ls"])
# os.system('java -jar BioTransformer3.0.jar  -k pred -b allHuman -ismi "CCCCCCCCC=CCCCCCCCCCCCC(=O)OC(COC(=O)CCCCCCCC=CCCCCCC)COP(=O)([O-])OCC[N+](C)(C)C" -ocsv phosphacho11.csv -s 1 -cm 1')

def biotransformer_para(smiles, chemical, s, cm):
    file_name = chemical + str(s) + str(cm)
    os.system(
        'java -jar BioTransformer3.0.jar  -k pred -b ecbased -ismi "' + smiles + '" -ocsv ' + file_name + '.csv -s ' + str(s) + ' -cm ' + str(cm))
    return file_name

def get_smiles_from_inchikey(inchikey):
    r = requests.get(
        f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/CanonicalSMILES,InChI/JSON')
    code = r.status_code
    print(code)
    if code != 200:
        csmiles = "Nah"
        cinchi = "Nah"
    else:
        r = r.json()
        csmiles = r['PropertyTable']['Properties'][0]['CanonicalSMILES']
        cinchi = r['PropertyTable']['Properties'][0]['InChI']
    return csmiles, cinchi


def inchi_list(file):
    inchi_l = []
    smile_l = []
    data = pd.read_csv(file)
    inchikeys = data["InChIKey"].values
    for i, inchikey in enumerate(inchikeys):
        s, i = get_smiles_from_inchikey(inchikey)
        inchi_l.append(i)
        smile_l.append(s)
        time.sleep(0.5)
        print(datetime.datetime.now().strftime("%H:%M:%S"))
    data["InChiKey2"] = inchikeys
    data["InChI2"] = inchi_l
    data["SMILES_new"] = smile_l
    return data

def hmdb_xml_match_merge(target_df, ref_df, filname):
    target_left = pd.read_csv(target_df, dtype=str)
    ref_right = pd.read_csv(ref_df, dtype=str)
    print(target_left.head(), ref_right.head())
    hmdb_matches = pd.merge(target_left, ref_right, how="left", on=["InChIKey"])
    hmdb_matches.to_csv(filname)
    print(hmdb_matches.head())
    return hmdb_matches

def biotransformer_loop(csv_file : str):
    os.chdir("biotransformer3.0jar")
    ref_file = "/Users/ciaraconway/Documents/HHEAR_BIO/hmdb_metabolites_final.csv"
    s = 2
    cm = 3
    csv_data = pd.read_csv(csv_file, dtype=str, encoding = "ISO-8859-1")
    smiles_list = csv_data["SMILES"].values
    for i, smiles in enumerate(smiles_list):
        chemical = csv_data.loc[csv_data.SMILES == smiles, "PREFERRED NAME"][i]
        chemical = chemical.replace(" ", "")
        filename = biotransformer_para(smiles, chemical, s, cm)
        filename_nist = "/Users/ciaraconway/biotransformer3.0jar/" + filename + ".csv"
        csv_nist = "/Users/ciaraconway/Documents/all_databases/Compounds/Pesticides/" + filename + "_Nist.csv"
        if os.path.exists(filename_nist):
            inchi_list(filename_nist).to_csv(csv_nist)
            filename_HMDB = "/Users/ciaraconway/Documents/all_databases/Compounds/Pesticides/" + chemical + str(s) + str(cm) + "_hmdb.csv"
            hmdb_xml_match_merge(csv_nist, ref_file, filename_HMDB)
        else:
            pass


###Call The Loop


def main():
    pass


if __name__ == '__main__':
    main()
