import glob
import os.path
import xlrd, openpyxl
import pandas as pd
from collections import defaultdict
import numpy as np


empty_df = pd.DataFrame()
csv_path = "/Users/ciaraconway/Documents/all_databases/NHANES_data/2017-2018/laboratory_data/"
path_msfinder = "/Users/ciaraconway/Documents/all_databases/NHANES_data/2017-2018/laboratory_data/*.xlsx"
count = 1
for i, lab_data in enumerate(glob.glob(path_msfinder)):
    if i == 0:
        df = pd.read_excel(lab_data, sheet_name='Sheet1')
        empty_df = df
        print(count)
        count += 1
    else:
        df_1 = pd.read_excel(lab_data, sheet_name='Sheet1')
        if len(df_1.index) > len(empty_df.index):
            empty_df = pd.merge(empty_df, df_1, on=["SEQN"], how="right")
            print(count)
            count += 1
        else:
            empty_df = pd.merge(empty_df, df_1, on=["SEQN"], how="left")
            print(count)
            count += 1

empty_df.to_csv(csv_path + "_all_NHAMES_" + ".csv")




