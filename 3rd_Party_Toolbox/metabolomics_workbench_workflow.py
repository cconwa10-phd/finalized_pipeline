import os
import subprocess
import datetime
import requests
import pandas as pd
import time

def parse_tsv_file(filetsv):
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


def hmdb_xml_match_merge_name(target_df, ref_df, filename): #GCMS Metabolomics Workbench HMDB Matches
    target_left = pd.read_csv(target_df, dtype=str, encoding = "ISO-8859-1")
    ref_right = pd.read_csv(ref_df, dtype=str, encoding = "ISO-8859-1")
    print(target_left.head(), ref_right.head())
    hmdb_matches = pd.merge(target_left, ref_right, how="left", on=["name"])
    hmdb_matches.to_csv(filename)
    print(hmdb_matches.head())
    return hmdb_matches

def metabolite_workbench_match(target_df, ref_df, filename): #GCMS Metabolomics Workbench Matches
    target_left = pd.read_csv(target_df, dtype=str, encoding="ISO-8859-1")
    ref_right = pd.read_csv(ref_df, dtype=str, encoding="ISO-8859-1")
    print(target_left.head(), ref_right.head())
    metabolite_match = pd.merge(target_left, ref_right, how="left", on=["pubchem_id"])
    metabolite_match.to_csv(filename)
    return metabolite_match


ref_file = "/Users/ciaraconway/Desktop/full_hmdb_xml_parsed_name.csv"
ref_file_m = "/Users/ciaraconway/Desktop/df_workbench.csv"

tsv = "/Users/ciaraconway/Desktop/st000058_new.tsv"
csv = "/Users/ciaraconway/Desktop/st000058_new.csv"
#parse_tsv_file(tsv).to_csv(csv)

filename_hmdb = "/Users/ciaraconway/Desktop/st000058_new_hmdb.csv"
#hmdb_xml_match_merge_name(csv, ref_file,filename_hmdb)

filename_m = "/Users/ciaraconway/Desktop/st000058_new_hmdb_MW.csv"
metabolite_workbench_match(filename_hmdb, ref_file_m, filename_m)


