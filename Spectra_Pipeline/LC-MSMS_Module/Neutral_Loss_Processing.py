import pandas as pd


def dbID_loop_msp(dbIDs, mspID):
    dbIDString = None
    for i, dbID in enumerate(dbIDs):
        if dbID == mspID:
            dbIDString = mspID
    return dbIDString

def neutralMassConverter(mspfile : str):
    count = 0
    neutralmsp = ""
    negfile = ""
    neutralmspFinal = ""
    precursormz = 0
    mz_values = []
    name = None

    smileCSV = "/Users/ciaraconway/Documents/MoNA - Clustering/pos_mode_smiles_cluster.csv"
    smileClusters = pd.read_csv(smileCSV)
    #dbIDs = smileClusters["DB"].values
    dbIDs = smileClusters.iloc[::100, :]
    dbIDs.to_csv("/Users/ciaraconway/Documents/MoNA - Clustering/pos_mode_smiles_cluster_100.csv")
    dbIDs = dbIDs["DB"].values


    for line in open(mspfile):
        line = line.strip()
        if len(line) > 0:
            if ":" in line:
                key, value = line.split(":", 1)
                if key.strip() == "PrecursorMZ":
                    value = value.strip()
                    precursormz += float(value.strip())
                    neutralmsp += line + "\n"
                elif key.strip() == "DB#":
                    name = dbID_loop_msp(dbIDs, value.strip())
                    neutralmsp += line + "\n"
                else:
                    neutralmsp += line + "\n"
            else:
                mz, intensity = map(float, line.split(' ', 1))
                mz_n = precursormz - mz
                if mz_n >= 0.001 and intensity >= 5:
                    mz_values.append(mz_n)
                    neutralmsp += str(mz_n) + " " + str(intensity) + "\n"
                #elif mz_n < 0 and intensity >= 5:
                    #neutralmsp += str(0) + " " + str(intensity) + "\n"
                #elif mz_n < -0.001:
                    #negfile += name + "\n" + str(mz_n) + "\t" + str(intensity) + "\n\n"
                    #print(name)
                #mz_values.append(mz_n)
                #neutralmsp += str(mz_n) + " " + str(intensity) + "\n"
        else:
            if name is not None:
                neutralmsp += "\n"
                neutralmspFinal += neutralmsp
                count += 1
                print(count)

            precursormz = 0
            mz_values = []
            name = None
            neutralmsp = ""
    return neutralmspFinal, negfile

def regMSPConvert(mspfile : str):
    count = 0
    msp = ""
    mspFinal = ""
    precursormz = 0
    mz_values = []
    intensities = []
    name = None

    smileCSV = "/Users/ciaraconway/Documents/MoNA - Clustering/pos_mode_smiles_cluster.csv"
    smileClusters = pd.read_csv(smileCSV)
    dbIDs = smileClusters.iloc[::100, :]
    dbIDs = dbIDs["DB"].values
    #dbIDs = dbIDs.iloc[::1000, :]

    for line in open(mspfile):
        line = line.strip()
        if len(line) > 0:
            if ":" in line:
                key, value = line.split(":", 1)
                if key.strip() == "PrecursorMZ":
                    value = value.strip()
                    precursormz += float(value.strip())
                    msp += line + "\n"
                elif key.strip() == "DB#":
                    name = dbID_loop_msp(dbIDs, value.strip())
                    msp += line + "\n"
                else:
                    msp += line + "\n"
            else:
                mz, intensity = map(float, line.split(' ', 1))
                if mz >= 0.001 and intensity >= 5:
                    mz_values.append(mz)
                    msp += str(mz) + " " + str(intensity) + "\n"
                #elif mz_n < 0 and intensity >= 5:
                    #neutralmsp += str(0) + " " + str(intensity) + "\n"
                #elif mz_n < -0.001:
                    #negfile += name + "\n" + str(mz_n) + "\t" + str(intensity) + "\n\n"
                    #print(name)
                #mz_values.append(mz_n)
                #neutralmsp += str(mz_n) + " " + str(intensity) + "\n"
        else:
            if name is not None:
                msp += "\n"
                mspFinal += msp
                count += 1
                print(count)

            precursormz = 0
            mz_values = []
            name = None
            msp = ""
    return mspFinal

def write_out_msp_new(mspFile, file):
    out_file = open(file, "w")
    out_file.write(mspFile)
    out_file.close()
    print("Outfile: ", file, " has been saved")



def main():
    print("Start")
    mspFile, negmspFile = neutralMassConverter("/Users/ciaraconway/Documents/MoNA - Clustering/MoNA-export-LipidBlast-09121.msp")
    print("file 1 done")
    file1 = "/Users/ciaraconway/Documents/MoNA - Clustering/MoNA-export-LipidBlast-09121-neutral2_neg_new_331.msp"
    write_out_msp_new(mspFile, file1)
    file2 = "/Users/ciaraconway/Documents/MoNA - Clustering/MoNA-export-LipidBlast-09121-reg_new_331.msp"
    regmspFile = regMSPConvert("/Users/ciaraconway/Documents/MoNA - Clustering/MoNA-export-LipidBlast-09121.msp")
    print("file 2 done")
    write_out_msp_new(regmspFile, file2)
    print("done")


if __name__ == '__main__':
    main()
