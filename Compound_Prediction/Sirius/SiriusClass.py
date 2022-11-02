import pandas as pd
import numpy as np
import os, shutil
import glob
from send2trash import send2trash
#os.environ(["PATH"])


def msp_to_msfile(mspFile, path):
    msFile = ""

    for line in open(mspFile):
        line = line.strip()
        if len(line) > 0:
            if ":" in line:
                key, value = line.split(":", 1)
                if key.strip() == "Name":
                    msFile += ">compound    " + value.strip() + "\n"
                elif key.strip() == "Precursor_type":
                    msFile += ">ionization  " + value.strip() + "\n"
                elif key.strip() == "PrecursorMZ":
                    msFile += ">parentmass  " + value.strip() + "\n"
                elif key.strip() == "Spectrum_type":
                    msFile += ">ms2" + "\n"
                # elif key.strip() == "Collision_energy":
                #     msFile += ">collision   " + value.strip() + "\n"
            else:
                msFile += line + "\n"
        else:
            msFile += "\n\n"
    outfile = open(path + "file1.ms", "w")
    outfile.write(msFile)
    outfile.close()
    return msFile

def run_sirius(msFile, output):
    #os.chdir("/Users/ciaraconway/")
    os.system("/Users/ciaraconway/sirius.app/Contents/MacOS/sirius --input " + msFile + " --output " + output + "out1" + " --ignore-formula formula --database=ALL -p default fingerprint structure compound-classes write-summaries --output " + output + "out2")
    #os.chdir("/Users/ciaraconway/new_complete_pipeline/")
def parse_compound_file(path):
    df = pd.read_table(path + "out2/compound_identifications.tsv", sep="\t")
    return df

def ms_sirius_parse(mspFile, path):
    msFile = ""
    df = pd.DataFrame()
    for line in open(mspFile):
        line = line.strip()
        if len(line) > 0:
            if ":" in line:
                key, value = line.split(":", 1)
                if key.strip() == "Name":
                    msFile += ">compound " + value.strip() + "\n"
                elif key.strip() == "Precursor_type":
                    msFile += ">ionization " + value.strip() + "\n"
                elif key.strip() == "PrecursorMZ":
                    msFile += ">parentmass " + value.strip() + "\n"
                elif key.strip() == "Spectrum_type":
                    msFile += ">ms2" + "\n"
                elif key.strip() == "Collision_energy":
                    msFile += ">collision " + value.strip() + "\n"
            else:
                msFile += line + "\n"
        else:
            outfile = open(path + "file1.ms", "w")
            outfile.write(msFile)
            outfile.close()
            run_sirius(path + "file1.ms", path)
            if os.path.exists(path + "out2/compound_identifications.tsv"):
                df = df.append(parse_compound_file(path))
                shutil.rmtree(path)
                os.mkdir(path)
            else:
                shutil.rmtree(path)
                os.mkdir(path)
            msFile = ""
    if df.empty:
        print("No Compound Predictions")
    else:
        return df

count = 0
dff = ms_sirius_parse("/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/cmpd_test_software.msp", "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/sirius_test/safe/")
dff.to_csv("/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/sirius_test/key_out.csv")
print(dff)
