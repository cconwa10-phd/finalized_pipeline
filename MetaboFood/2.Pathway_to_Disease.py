import collections

import pandas as pd
import numpy as np

#Maps the KEGG/Reactome IDS to the Disease from CTD - Need to have CTD_diseases_pathways.csv
def comparing_files(kayFile, diseaseFile):
    disDict = collections.defaultdict(list)
    kay = pd.read_csv(kayFile)
    disease = pd.read_csv(diseaseFile)
    count = 0
    kayPathIDs = kay["IDS"].values
    kayPathIDs = kayPathIDs.unique()
    for i, kayPathID in enumerate(kayPathIDs):
        disPathIDs = disease["PathwayIDREAL"].values
        for j, disPathID in enumerate(disPathIDs):
            if disPathID.strip() == kayPathID.strip():
                disDict["Kay ID"].append(kayPathID)
                disDict["Disease ID"].append(disPathID)
                disDict["Disease Name"].append(disease.loc[disease.PathwayIDREAL == disPathIDs, "DiseaseName"][j])
                disDict["MESH"].append(disease.loc[disease.PathwayIDREAL == disPathIDs, "DiseaseID"][j])
                disDict["Pathway ID"].append(disease.loc[disease.PathwayIDREAL == disPathIDs, "PathwayName"][j])
                disDict["Gene ID"].append(disease.loc[disease.PathwayIDREAL == disPathIDs, "InferenceGeneSymbol"][j])
                print(count)
                count += 1
    print(disDict)
    df = pd.DataFrame.from_dict(disDict)
    df.to_csv("/Users/cconwa10/Desktop/disease_kegg_path_match.csv")
    print(df.head(10))
    return df


def main():
    pass


if __name__ == '__main__':
    main()
