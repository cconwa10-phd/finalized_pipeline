from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

#!pip install scipy
#!pip install matplotlib
#!pip install numpy

import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
import matplotlib.pyplot as plt
import matplotlib.style as style
import numpy as np
from rdkit import DataStructs
from rdkit.Chem import Draw
#from rdkit.Chem.Draw import IPythonConsole

import sys
import os

def batch_smile_validation(file):
    smilesdf = pd.read_csv(file, dtype=str, encoding = "ISO-8859-1")
    #smilesdf = smilesdf[smilesdf["SMILES"] != "Nah"]
    smiles = smilesdf["CanonicalSMILES"].values
    valid_list = []
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        print(mol)
        valid_list.append(mol)
    smilesdf["Check"] = valid_list
    smilesdf.to_csv(file)
    return smilesdf

def single_smile_validation(smile):
    mol = Chem.MolFromSmiles(smile)
    pass

def structure_similarity_search(tSmile, refSmiles, tolerance, file, pic_label, rfps):
    new_list = []
    new_scores = []
    mols = []
    # refSmiles = refSmiles.drop_duplicates(subset=["InChIKey"])
    # refSmiles.reset_index(drop=True)
    rSmiles = refSmiles["SMILES"].values
    tMol = Chem.MolFromSmiles(tSmile)
    tFPS = AllChem.GetMorganFingerprintAsBitVect(tMol,2,1024)

    for i, rfp in enumerate(rfps):

        fp1 = tFPS
        fp2 = rfp

        similarity = DataStructs.TanimotoSimilarity(fp1,fp2)
        if similarity >= tolerance:
            row = refSmiles.iloc[[i]].values.flatten().tolist()
            new_list.append(row)
            new_scores.append(similarity)
            rMol = Chem.MolFromSmiles(rSmiles[i])
            mols.append(rMol)
    new_dataframe = pd.DataFrame(data=new_list, columns=refSmiles.columns)
    new_dataframe["mols"] = mols
    new_dataframe["Score"] = new_scores
    new_dataframe.sort_values(by=["Score"])
    new_dataframe.to_csv(file)
    if os.path.exists(file):
        generating_sds_files(pic_label, new_dataframe)
    else:
        print("None")


def generating_sds_files(tLabel, nrefDataframe):
    file_label = tLabel
    molL = nrefDataframe['mols'].tolist()
    labelL = nrefDataframe['name'].tolist()
    scoreL = nrefDataframe['Score'].tolist()
    labelF = list(zip(labelL,scoreL))
    list_length = len(molL)
    counter = 1
    tempmolL = []
    for i, mol in enumerate(molL):
        if len(tempmolL) < 50 and list_length != 0:
            tempmolL.append(mol)
            list_length = list_length - 1
        elif len(tempmolL) == 50:
            mol_show = Draw.MolsToGridImage(molL, molsPerRow=3, subImgSize=(500, 500), legends=[str(i) for i in labelF],
                                            returnPNG=False)
            file_label = file_label + str(counter)
            genpath = "/Users/ciaraconway/Desktop/fig_HHEAR/"
            filepath = os.path.join(genpath, 'generated_cluster_{}.png'.format(file_label))
            mol_show.save(filepath)
            tempmolL = []
            counter += 1
    mol_show = Draw.MolsToGridImage(molL, molsPerRow=3, subImgSize=(500, 500), legends=[str(i) for i in labelF],
                                    returnPNG=False)
    file_label = file_label + str(counter)
    genpath = "/Users/ciaraconway/Desktop/fig_HHEAR/"
    filepath = os.path.join(genpath, 'generated_cluster_{}.png'.format(file_label))
    mol_show.save(filepath)

def looping_target_smiles(tDataFrame, refDataFrame, tolerance):
    file_path = "/Users/ciaraconway/Desktop/HHEAR_matches/"
    tSmiles = pd.read_csv(tDataFrame, dtype=str, encoding="ISO-8859-1")
    cSmiles = tSmiles["CanonicalSMILES"].values
    tLabel = tSmiles["Full Analyte Name"].values
    refSmiles = pd.read_csv(refDataFrame, dtype=str, encoding="ISO-8859-1")
    rSmiles = refSmiles["SMILES"].values
    rMol = [Chem.MolFromSmiles(smile) for smile in rSmiles]
    rfps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024) for mol in rMol]
    for i, smile in enumerate(cSmiles):
        name = 'matches_{}.csv'.format(tLabel[i])
        file = file_path + name
        structure_similarity_search(smile,refSmiles,tolerance, file, tLabel[i], rfps)





refDataFrame = "/Users/ciaraconway/Documents/all_databases/all_ref_smiles_valid.csv"
tDataFrame = "/Users/ciaraconway/Documents/BioTransTest/Databases/all_LCMSMS.csv"
#batch_smile_validation(tDataFrame)

looping_target_smiles(tDataFrame, refDataFrame, 0.95)

