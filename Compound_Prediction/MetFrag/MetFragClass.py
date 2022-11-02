import numpy as np
import pandas as pd
import os


def write_parameter_file(peak_path, ion_mass, result_path, pos_mode):
    para_file = ""
    para_file += "PeakListPath = " + peak_path + "\n" + \
                 "IonizedPrecursorMass = " + ion_mass + "\n" + \
                 "ResultsPath = " + result_path + "\n" + \
                 "IsPositiveIonMode = " + pos_mode + "\n"

    para_file += '''MetFragDatabaseType = PubChem
DatabaseSearchRelativeMassDeviation = 5
FragmentPeakMatchAbsoluteMassDeviation = 0.008
FragmentPeakMatchRelativeMassDeviation = 5
PrecursorIonMode = 1
MetFragScoreTypes = FragmenterScore
MetFragScoreWeights = 1.0
MetFragCandidateWriter = CSV
SampleName = results
MaximumTreeDepth = 2
MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter
MetFragPostProcessingCandidateFilter = InChIKeyFilter'''
    out_file = open("parameter_file.txt", "w")
    out_file.write(para_file)
    out_file.close()
    return para_file


def parse_peaks(msp_file):
    count = 0
    cwd = "/Users/ciaraconway/new_complete_pipeline/Compound_Prediction/MetFrag/"
    df = pd.DataFrame()
    is_present = ["FALSE", "FALSE", "FALSE"]
    params = []
    peaks = ""
    for line in open(msp_file):
        line = line.strip()
        if len(line) > 0:
            if ":" in line:
                key, value = line.split(":", 1)
                if key.strip() == "Name":
                    params.append(value.strip())
                    is_present[0] = "TRUE"
                elif key.strip() == "PrecursorMZ":
                    params.append(value.strip())
                    is_present[1] = "TRUE"
                elif key.strip() == "Ion_mode":
                    if value.strip() == "P":
                        params.append("True")
                        is_present[2] = "TRUE"
                    else:
                        params.append("False")
                        is_present[2] = "TRUE"
            else:
                mz, intensity = line.split("\t", 1)
                peaks += mz.strip() + "\t" + intensity.strip() + "\n"
        elif "FALSE" in is_present:
            pass
            params = []
            peaks = ""
            count += 1
            print(count)
            print(is_present)
            is_present = ["FALSE", "FALSE", "FALSE"]
        else:
            out_file = open("peaks.txt", "w")
            out_file.write(peaks)
            out_file.close()
            write_parameter_file("peaks.txt", str(params[1]), cwd, params[2])
            os.system("java -jar /Users/ciaraconway/MetFragCommandLine-2.5.0.jar parameter_file.txt")
            r_df = pd.read_csv("results.csv")
            r_df = r_df.head(20)
            print(r_df.shape[0])
            r_df["Name"] = [params[0]]*r_df.shape[0]
            df = df.append(r_df)
            count += 1
            print(count)
            print(is_present)
            params = []
            peaks = ""
            is_present = ["FALSE", "FALSE", "FALSE"]
    return df


file = "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/cmpd_test_software.msp"
parse_peaks(file).to_csv("final_results.csv")
