import datetime
import requests
import pandas as pd
import time



def inch_loop(data_inchks, msp_inchk):
    found = False
    for i, inchik in enumerate(data_inchks):
        if inchik.strip() == msp_inchk:
            found = True
    return found

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def assign_spectra(filename, msp):
    count = 0
    mspFile = ""
    mzval = []
    intval = []
    meta_data = {}
    data = pd.read_csv(filename, encoding = "ISO-8859-1")
    inchiks = data["InChIKey_use"].values

    count = 0
    for line in open(msp):
        count += 1
        line = line.strip()
        if len(line) > 0:
            if ":" in line:
                key, value = line.split(":", 1)
                if key.strip() == "Name":
                    value = value.strip() + ":::" + str(count)
                    meta_data[key] = value
                else:
                    meta_data[key] = value
            else:
                mz, intensity = line.split(' ', 1)
                mzval.append(float(mz.strip()))
                if isfloat(intensity):
                    intval.append(float(intensity.strip()))
                else:
                    inten, extra = intensity.split(" ", 1)
                    intval.append(float(inten.strip()))
        else:
            if meta_data.get("InChIKey") is not None:
                found = inch_loop(inchiks, meta_data.get("InChIKey").strip())
                if found:
                    count += 1
                    for kMeta, vMeta in meta_data.items():
                        mspFile += str(kMeta) + ": " + str(vMeta) + "\n"
                    for kMZ, vInt in zip(mzval, intval):
                        mspFile += str(kMZ) + "\t" + str(vInt) + "\n"
                    mspFile += "\n"
                    print(meta_data)
                else:
                    pass
            else:
                pass
            mzval = []
            intval = []
            meta_data = {}
    #print(mspFile)
    print(count)
    return mspFile

def write_out_msp_USDA(mspFile, file):
    out_file = open(file, "w")
    out_file.write(mspFile)
    out_file.close()
    print("Outfile: ", file, " has been saved")



def main():
    msp_file_search = "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/lcmsmspn_MoNA.msp"
    csv_match = "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/cmpd_test_software.csv"
    new_msp = assign_spectra(csv_match, msp_file_search)
    out = "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/cmpd_test_software.msp"
    write_out_msp_USDA(new_msp, out)

if __name__ == '__main__':
    main()
