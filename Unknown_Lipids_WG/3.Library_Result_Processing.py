import glob
import os.path
import datetime

import pandas as pd
from collections import defaultdict
import numpy as np
import requests
import time
import urllib.request


###identity search
#MSPepSearch64.exe zD /Z 0.01 /M 0.01 /MatchPolarity /LIB "C:\NIST20\MSSEARCH\hr_msms_nist" /LIB "C:\NIST20\MSSEARCH\lr_msms_nist" /LIB "C:\NIST20\MSSEARCH\apci_msms_nist" /LIB "C:\NIST20\MSSEARCH\biopep_msms_nist"  /LIB "C:\NIST20\MSSEARCH\MoNA-export-LipidBlast” /MinMF 300 /MinInt 5 /HITS 10 /OutChemForm /OutCASrn /OutMW /OutIK /OutMaxScore /INP "C:\Users\cconwa10\Desktop\MSMS_UCD_POS\NEG_MSMS.msp" /OUTTAB "C:\Users\cconwa10\Desktop\MSMS_UCD_POS\msms_UCD_pos_identity.tsv"

###hybrid search
#MSPepSearch64.exe yD /Z 0.01 /M 0.01 /MatchPolarity /LIB "C:\NIST20\MSSEARCH\hr_msms_nist" /LIB "C:\NIST20\MSSEARCH\lr_msms_nist" /LIB "C:\NIST20\MSSEARCH\apci_msms_nist" /LIB "C:\NIST20\MSSEARCH\biopep_msms_nist"  /LIB "C:\NIST20\MSSEARCH\MoNA-export-LipidBlast” /MinMF 300 /MinInt 5 /HITS 10 /OutChemForm /OutCASrn /OutMW /OutIK /OutMaxScore /INP "C:\Users\cconwa10\Desktop\MSMS_UCD_POS\NEG_MSMS.msp" /OUTTAB "C:\Users\cconwa10\Desktop\MSMS_UCD_POS\msms_UCD_pos_hybrid.tsv"

####will do this twice - one for hybrid, one for exact match
def parse_tsv_NIST(filetsv):
    nistResults = open(filetsv, "r")
    processedResults = []
    processedColumns = []
    lineCount = 0

    inchiL = []
    smilesL = []

    for line in nistResults:
        parsedInfo = line.split("\t")
        if (lineCount == 3) & (len(parsedInfo) > 5):
            print(parsedInfo)
            headers = parsedInfo
            processedColumns = headers
        elif (lineCount >= 4) & (len(parsedInfo) > 5):
            print(parsedInfo)
            processedResults.append(parsedInfo)
        lineCount += 1
    df = pd.DataFrame(processedResults, columns=processedColumns)
    print(df.head(10))

    inchikeys = df["InChIKey"].values
    for inchikey in inchikeys:
        inchiN = inchikey_to_inchi(inchikey)
        smilesN = inchi_to_smiles(inchiN)
        inchiL.append(inchiN)
        smilesL.append(smilesN)
        #df["InChI_U"] = inchiU
    df["InChI"] = inchiL
    df["SMILES"] = smilesL
    print(df.head(10))
    return df

def parse_tsv_NIST_only(filetsv):
    nistResults = open(filetsv, "r")
    processedResults = []
    processedColumns = []
    lineCount = 0

    inchiL = []
    smilesL = []

    for line in nistResults:
        parsedInfo = line.split("\t")
        if (lineCount == 3) & (len(parsedInfo) > 5):
            print(parsedInfo)
            headers = parsedInfo
            processedColumns = headers
        elif (lineCount >= 4) & (len(parsedInfo) > 5):
            print(parsedInfo)
            processedResults.append(parsedInfo)
        lineCount += 1
    df = pd.DataFrame(processedResults, columns=processedColumns)
    return df

###will parse msfinder file
def parse_MSFINDER(filetxt, samplename):
    count = 0
    parsedInfo = []
    infoHold = []
    columns = ["Name", "InchiKey", "Smiles", "TotalScore", "SampleName"]
    msfinderResults = open(filetxt, "r")
    for line in msfinderResults:
        if len(line) > 0:
            if ":" in line:
                key,value = line.split(":", 1)
                if key.strip() == "NAME":
                    infoHold.append(value.strip())
                elif key.strip() == "INCHIKEY":
                    infoHold.append(value.strip())
                elif key.strip() == "SMILES":
                    infoHold.append(value.strip())
                elif key.strip() == "TotalScore":
                    infoHold.append(value.strip())
                elif key.strip() == "CcsSimilarityScore": ##this is used as a marker to append the name
                    infoHold.append(samplename)
                    parsedInfo.append(infoHold)
                    infoHold = []
            else:
                pass
        else:
            pass
    df = pd.DataFrame(parsedInfo, columns=columns)
    print(df.head(10))
    return df
###will search ChemSpider for inchi, will provide csv file
def inchikey_to_inchi(inchik):
    host = "http://www.chemspider.com"
    getstring = "/InChI.asmx/InChIKeyToInChI?inchi_key="
    inchikey = inchik

    r = requests.get('{}{}{}'.format(host, getstring, inchikey))
    if r.ok:
        content = r.text
        res = str(r.text.replace('<?xml version="1.0" encoding="utf-8"?>\r\n<string xmlns="http://www.chemspider.com/">', '').replace('</string>', '').strip())
        print(res)
    else:
        res = "Nah"
    return res

def inchi_to_smiles(inchi):
    host = "http://www.chemspider.com"
    getstring = "/InChI.asmx/InChIToSMILES?inchi="
    inch = str(inchi.strip())

    r = requests.get('{}{}{}'.format(host, getstring, inch))
    if r.ok:
        res = str(r.text.replace('<?xml version="1.0" encoding="utf-8"?>\r\n<string xmlns="http://www.chemspider.com/">', '').replace('</string>', '').strip())
        print(res)
    else:
        res = "Nah"
    return res


def get_smiles_from_inchikey(inchikey):
    r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/'
                     f'compound/inchikey/{inchikey}/property/'
                     f'CanonicalSMILES,InChI/JSON')
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
    data["SMILES2"] = smile_l
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


def main():
    # file = "/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Files/all_unk/newfile_WG_unknown_all_n_hybrid_700_20.tsv"
    # result = "/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Files/all_unk/newfile_WG_unknown_all_n_hybrid_700_20.csv"
    # parse_tsv_NIST(file).to_csv(result)
    file1 = "/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Files/neg_/unk_neg_spectra_hybrid_700_20.tsv"
    result1 = "/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Files/neg_/unk_neg_spectra_hybrid_700_20.csv"
    result2 = "/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Files/neg_/unk_neg_spectra_hybrid_700_20_p.csv"
    #parse_tsv_NIST(file1).to_csv(result1)
    parse_tsv_NIST_only(file1).to_csv(result1)
    inchi_list(result1).to_csv(result2)
    # empty_df = pd.DataFrame(columns=["Name", "InchiKey", "Smiles", "TotalScore", "SampleName"])
    # csv_path = "/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSFINDER/MSFINDER_CSV/"
    # path_msfinder = "/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSFINDER/MSFINDER_NEW_N/*"
    # count = 0
    # for fname in glob.glob(path_msfinder):
    #     path_text = fname+"/*.sfd"
    #     for sfdname in glob.glob(path_text):
    #         filesize = os.path.getsize(sfdname)
    #         if filesize == 0:
    #             pass
    #         else:
    #             name = fname.split("/")[7]
    #             #filled_df = parse_MSFINDER(sfdname, name)
    #             empty_df = empty_df.append(parse_MSFINDER(sfdname, name))
    # count += 1
    # empty_df.to_csv(csv_path + "_NEW_N_" + str(count) + ".csv")


if __name__ == '__main__':
    main()
