import pandas as pd
import numpy as np


def pos_neg_compare(pos, neg, file):
    posdf = pd.read_csv(pos, usecols=[0,1,2], dtype = str)
    negdf = pd.read_csv(neg, dtype = str)
    posdf = posdf.drop_duplicates("Alignment_ID")
    unk_matches = pd.merge(posdf, negdf, how="left", on=["Alignment_ID"])
    unk_matches.to_csv(file)


path = "/Users/ciaraconway/Documents/MSMS_lipids_allLabs/"
pos = path + "all_unk_classfied_version_only.csv"
neg = path + "unk_neg_csv.csv"
file = path + "unk_neg_pos_matches.csv"

pos_neg_compare(pos, neg, file)
