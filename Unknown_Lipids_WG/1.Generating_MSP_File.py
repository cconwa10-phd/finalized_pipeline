import os, sys
import glob
import csv
import pandas as pd
import numpy as np

def calc_exactmass_new(precursor, adduct):
    #proton = 0.000549
    if adduct.strip() == "[M+2H]2+":
        exactMass = precursor/2 - 1.007276
    elif adduct.strip() == "[M+H]+":
        exactMass = precursor - 1.007276
    elif adduct.strip() == "[M+Na]+":
        exactMass = precursor - 22.989218
    elif adduct.strip() == "[M+NH4]+":
        exactMass = precursor - 18.033823
    elif adduct.strip() == "[M+H-H2O]+":
        exactMass = precursor - 1.007276 + 18.010565
    elif adduct.strip() == "[M+H]2+":
        exactMass = (precursor - 1.007276)/2
    elif adduct.strip() == "[M+Na]2+":
        exactMass = (precursor - 22.989218)/2
    else:
        exactMass = None
    return exactMass

def read_unknowns_new(filenameCSV, filenameMSP):
    ###Tolerances
    mzTolerance = 0.008
    rtTolerance = 0.1

    ###Reading in Unknown VIP CSV
    top50Unknowns = pd.read_csv(filenameCSV)
    averageMZs = top50Unknowns["Average_Mz"].values

    ###Creating MSP
    mspFile = ""
    metaData = {}
    mzValues = []
    intValues = []

    for i, averageMZ in enumerate(averageMZs):
        count = 0
        alignmentId = top50Unknowns.loc[top50Unknowns.Average_Mz == averageMZ, "Alignment_ID"][i]
        rt = float(top50Unknowns.loc[top50Unknowns.Average_Mz == averageMZ, "Average_Rt"][i])
        for line in open(filenameMSP):
            line = line.strip()
            if len(line) > 0:
                if ":" in line:
                    key, value = line.split(":", 1)
                    if key.strip() == "NAME":
                        metaData["Name"] = value.strip() + ":::" + str(alignmentId)
                    elif key.strip() == "PRECURSORMZ":
                        metaData["PrecursorMZ"] = float(value.strip())
                    elif key.strip() == "PRECURSORTYPE":
                        pass
                    elif key.strip() == "Comment":
                        pass
                    else:
                        metaData[key] = value.strip()
                else:
                    mz, intensity = line.split("\t", 1)
                    mzValues.append(float(mz.strip()))
                    intValues.append(float(intensity.strip().split("\n")[0]))
            else:
                if metaData.get("Num Peaks") != str(0):
                    precur = float(metaData.get("PrecursorMZ"))
                    if abs(precur - float(averageMZ)) <= mzTolerance:
                        rtm = float(metaData.get("RETENTIONTIME"))
                        if abs(rtm - rt) <= rtTolerance:
                            count += 1
                            print("Found: Count = " + str(count) + " Alignment ID " + str(alignmentId))
                            metaData["IONMODE"] = "P"
                            metaData["Precursortype"] = top50Unknowns.loc[top50Unknowns.Average_Mz == averageMZ, "Adduct_type"][i]
                            metaData["Name"] = str(metaData["Name"]) + "^" + str(count)
                            print(metaData["Name"])
                            precurtype = metaData.get("Precursortype")
                            metaData["MW"] = calc_exactmass_new(precur, precurtype)
                            del metaData["Num Peaks"]
                            metaData["Num Peaks"] = len(mzValues)
                            for kMeta, vMeta in metaData.items():
                                mspFile += str(kMeta) + ": " + str(vMeta) + "\n"
                            for kMZ, vInt in zip(mzValues, intValues):
                                mspFile += str(kMZ) + "\t" + str(vInt) + "\n"
                            mspFile += "\n"
                    mzValues = []
                    intValues = []
                    metaData = {}
    print(mspFile)
    return mspFile

def write_out_msp_new(mspFile, file):
    out_file = open(file, "w")
    out_file.write(mspFile)
    out_file.close()
    print("Outfile: ", file, " has been saved")

def main():
    rawFile = '/Users/ciaraconway/Documents/MSMS_lipids_allLabs/Data_Dump_Tong/UnkLip_pos_MS2spectra.msp'
    top50 = '/Users/ciaraconway/Documents/MSMS_lipids_allLabs/Data_Dump_Tong/knowns_pos.csv'
    fileNew = '/Users/ciaraconway/Documents/MSMS_lipids_allLabs/Data_Dump_Tong/knowns_pos_found.msp'
    mspFile = read_unknowns_new(top50, rawFile)
    write_out_msp_new(mspFile, fileNew)




if __name__ == '__main__':
    main()
