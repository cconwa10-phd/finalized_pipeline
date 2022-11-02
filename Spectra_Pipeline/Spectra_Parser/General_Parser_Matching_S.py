import pandas as pd
import numpy as np

def parse_whole_msp_to_csv(msp):
    pass # we aren't doing this right now with your stuff

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
    mspFile = ""
    mzval = []
    intval = []
    meta_data = {}
    data = pd.read_csv(filename, encoding = "ISO-8859-1")
    inchiks = data["InChIKey_use"].values

    for line in open(msp):
        line = line.strip()
        if len(line) > 0:
            if ":" in line:
                key, value = line.split(":", 1)
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
    print(mspFile)
    return mspFile

def write_out_msp_(mspFile, file):
    out_file = open(file, "w")
    out_file.write(mspFile)
    out_file.close()
    print("Outfile: ", file, " has been saved")


def main():
    pass ### run the functions you want to run here
    msp_new = assign_spectra("LCMSMS.msp")
    write_out_msp_(msp_new, "new_file.msp")

if __name__ == '__main__':
        main()
