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

def smile_validation(file):
    smilesdf = pd.read_csv(file, dtype=str, encoding = "ISO-8859-1")
    #smilesdf = smilesdf[smilesdf["SMILES"] != "Nah"]
    smiles = smilesdf["SMILES"].values
    valid_list = []
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        print(mol)
        valid_list.append(mol)
    smilesdf["Check"] = valid_list
    smilesdf.to_csv("/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/HHEAR_Results/_run_validate.csv")
    return smilesdf

def ClusterFpsS(fps, tol):
    size = len(fps)
    if size > 1:
        #print(size)
        similarityM = np.zeros((size, size))

        for i in range(size):
            for j in range(size):
                fp1 = fps[i]
                fp2 = fps[j]
                similarity = DataStructs.TanimotoSimilarity(fp1,fp2)
                similarityM[i, j] = similarity
        distanceM = 1 - similarityM
        ndistanceM = ssd.squareform(distanceM, force = 'tovector', checks = True)
        Z = sch.linkage(ndistanceM, 'complete')
        #dn = sch.dendrogram(Z, leaf_rotation=270, leaf_font_size=3, labels = labels, color_threshold=0.15)

        hc = sch.fcluster(Z, t = tol, criterion = 'distance', depth=2, R=None, monocrit=None)

        cs = hc.tolist()
    else:
        cs = 1
    #print(len(cs))
    return cs

def alignment_grouping(file, tol):
    count_lists = []
    df = pd.read_csv(file, index_col=0, header=0)
    ids = df['Version'].tolist()
    unique_ids = set(ids)
    groups = [df[df.Version == uid] for uid in unique_ids]
    for i, group in enumerate(groups):
        smilesdf = group
        file_label = str(smilesdf["Version"].values[0])
        smiles = smilesdf["SMILES"].values
        labels = smilesdf["Peptide"].values
        mols = [Chem.MolFromSmiles(smile) for smile in smiles]
        fps = [AllChem.GetMorganFingerprintAsBitVect(mol,2,1024) for mol in mols]
        clustersS = ClusterFpsS(fps, tol)
        smilesdf["mols"] = mols
        smilesdf["clusters"] = clustersS
        smilesdf = smilesdf.sort_values(by=["clusters"])
        smilesdfc = smilesdf["clusters"].value_counts().sort_values(axis=0, ascending=False)
        print(smilesdfc.head(10))
        sumclust1 = smilesdfc.iloc[0:1].sum()/len(smilesdf["clusters"].values)
        sumclust2 = smilesdfc.iloc[0:2].sum()/len(smilesdf["clusters"].values)
        sumclust = smilesdfc.iloc[0:3].sum()/len(smilesdf["clusters"].values)
        count_lists.append((file_label, len(smilesdf["clusters"].values), sumclust1, sumclust2, sumclust))
        molL = smilesdf['mols'].tolist()
        labelL = smilesdf['Type'].tolist()
        clusterL = smilesdf['clusters'].tolist()
        #massL = smilesdf['Mass'].tolist()
        scoreL = smilesdf['Score'].tolist()
        # labelF = list(zip(labelL,scoreL,clusterL))
        # list_length = len(molL)
        # counter = 1
        # tempmolL = []
        # for i, mol in enumerate(molL):
        #     if len(tempmolL) < 50 and list_length != 0:
        #         tempmolL.append(mol)
        #         list_length = list_length - 1
        #     elif len(tempmolL) == 50:
        #         mol_show = Draw.MolsToGridImage(molL, molsPerRow=3, subImgSize = (500, 500), legends = [str(i) for i in labelF], returnPNG = True)
        #         file_label = file_label + str(counter)
        #         genpath = "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/HHEAR_Results/figs/"
        #         filepath = os.path.join(genpath, 'generated_cluster_{}.png'.format(file_label))
        #         mol_show.save(filepath)
        #         tempmolL = []
        #         counter += 1
        # mol_show = Draw.MolsToGridImage(molL, molsPerRow=3, subImgSize = (500, 500), legends = [str(i) for i in labelF], returnPNG = False)
        # file_label = file_label + str(counter)
        # genpath = "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/HHEAR_Results/figs/"
        # filepath = os.path.join(genpath, 'generated_cluster_{}.png'.format(file_label))
        # mol_show.save(filepath)
        # tempmolL = []
    df_counts = pd.DataFrame(count_lists, columns=["Name", "Length", "Ratio_1", "Ratio_2", "Ratio_3"])
    dfNew = pd.concat(groups)
    #df_counts.to_csv("/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/HHEAR_Results/all_results_count.csv")
    #dfNew.to_csv("/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/HHEAR_Results/all_results_clusters.csv")
    return dfNew
def rdkit_fig(file):
    pass


def main():
    file = "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/HHEAR_Results/all_results_less.csv"
    #smile_validation(file)
    file2 = "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/HHEAR_Results/_run_validate.csv"
    empty_data = pd.DataFrame()
    tols = [0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
    for tol in tols:
        df = alignment_grouping(file2, tol)
        cluster = df["clusters"].tolist()
        id = df['Version'].tolist()
        typ = df["Type"].tolist()
        empty_data[str(tol) + "version"] = id
        empty_data[str(tol) + "type"] = typ
        empty_data[str(tol) + "clusters"] = cluster
        print(empty_data)
    empty_data.to_csv("/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/HHEAR_Results/fm_hhear_data.csv")
if __name__ == '__main__':
    main()
