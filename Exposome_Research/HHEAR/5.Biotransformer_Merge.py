import pandas as pd
import numpy as np
import glob
import os

def comb_spreadsheets(path, files):
    joined_files = os.path.join(path, files)
    joined_list = glob.glob(joined_files)

    df = pd.concat(map(pd.read_csv, joined_list), ignore_index=True)
    return df



path = "/Users/ciaraconway/Documents/all_databases/Compounds/Pesticides/biotransformations"
files = "*hmdb.csv"

comb_spreadsheets(path, files).to_csv("/Users/ciaraconway/Documents/all_databases/Compounds/Pesticides/all_biotrans.csv")
