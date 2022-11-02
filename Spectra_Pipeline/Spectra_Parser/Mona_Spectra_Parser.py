import collections

import pandas as pd
import numpy as np


def lipidblast_parser(file):
    count = 0
    df = pd.DataFrame()
    meta_lipids = {}
    for line in open(file):
        if len(line) > 0:
            if ":" in line:
                key, value = line.split(":", 1)
                if key.strip() == "Comments":
                    value = value.split("\"")
                    for val in value[1::2]:
                        k,v = val.split("=", 1)
                        meta_lipids[k.strip()] = v.strip()
                elif key.strip() == "Num Peaks":
                    meta_lipids[key] = str(value.strip())
                    df = df.append(meta_lipids, ignore_index=True)
                    count += 1
                    print(count)
                    meta_lipids = {}
                else:
                    meta_lipids[key.strip()] = str(value.strip())
            else:
                pass
        else:
            pass
    df = df.drop_duplicates(subset="InChIKey").set_index("Name")
    print(df.head(10))
    return df



#full_parse_lipidblast("/Users/cconwa10/Documents/MoNA-export-LipidBlast.msp").to_csv("/Users/cconwa10/Documents/MoNA-export-LipidBlast-smiles-full.csv")
#parse_mz_intensities("/Users/ciaraconway/Documents/MoNA - Clustering/MoNA-export-LipidBlast-09121-neutral2_neg0.msp")


