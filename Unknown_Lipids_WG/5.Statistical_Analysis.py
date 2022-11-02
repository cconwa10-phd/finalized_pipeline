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


def alignment_grouping_new(file):
    df_list = []
    df = pd.read_csv(file, index_col=0, header=0)
    ids = df["Version"].tolist()
    unique_ids = set(ids)
    groups = [df[df.Version == uid] for uid in unique_ids]
    for i, group in enumerate(groups):
        smilesdf = group
        #file_label = str(smilesdf["Version"].values[0])
        #cluster = smilesdf["clusters"]
        smilesdf = smilesdf.sort_values(by=["clusters"]).iloc[:, 28:49]
        smilesdfd = smilesdf.drop_duplicates(subset="clusters")
        smilesdfc = smilesdf["clusters"].value_counts().sort_values(axis=0, ascending=False)
        smilesdfc = smilesdfc.to_frame(name="count")
        smilesdfc.reset_index(inplace=True)
        smilesdfc.rename(columns = {"index":"clusters"}, inplace=True)
        #smilesdfc = smilesdfc.columns = ["clusters", "count_clusters"]
        if len(smilesdfc) >= 3: #change between top 1 and top 3
            smilesdfc = smilesdfc.iloc[0:3, :]
        else:
            pass
        cluster_matches = pd.merge(smilesdfc, smilesdfd, how="left", on=["clusters"])
        cluster_matches.head(10)
        df_list.append(cluster_matches)
    results = pd.concat(df_list)
    results.to_csv("/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Files/ALL_EVERYTHING/known_top3n.csv")


def match_tiss_alignment(tiss):
    df_t = pd.read_csv(tiss)
    col = df_t.columns.values.tolist()
    col = col[:-1]
    df_t = (df_t.set_index(col).apply(lambda x: x.str.split(' ').explode()).reset_index())

    return df_t





def main():
    # file = "/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Files/Priority_Unkns/one_tissue_priority_results.csv"
    # smile_validation(file)
    file2 = "/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Files/ALL_EVERYTHING/known_n.csv"
    alignment_grouping_new(file2)

    #file3 = "/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Files/ALL_EVERYTHING/pos_top3.csv"
    #match_tiss_alignment(file3).to_csv("/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Files/ALL_EVERYTHING/pos_top3_tiss.csv")

if __name__ == '__main__':
    main()
