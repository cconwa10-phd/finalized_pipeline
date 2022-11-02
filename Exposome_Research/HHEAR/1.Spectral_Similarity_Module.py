import numpy as np
import pandas as pd
import pickle
from spectrum import Spectrum
from typing import List
import scipy.sparse as ss
import multiprocessing as mp


def remove_peaks(mz_values, intensities):
    mz_values_keep = []
    intensities_keep = []
    for i in range(len(mz_values)):
        mz = mz_values[i]
        intensity = intensities[i]
        if intensity >= 5:
            mz_values_keep.append(mz)
            intensities_keep.append(intensity)
    return mz_values_keep, intensities_keep


def remove_peaks_P_N(mz_values, intensities, exclude_mz):
    mz_tolerance = 0.01
    mz_values_keep = []
    intensities_keep = []
    for i in range(len(mz_values)):
        mz = mz_values[i]
        intensity = intensities[i]
        if abs(mz - exclude_mz) >= mz_tolerance:
            mz_values_keep.append(mz)
            intensities_keep.append(intensity)
    return mz_values_keep, intensities_keep


def read_file(filename):
    index = 0
    spectra_p = []
    spectra_n = []
    name = None
    mz_values = []
    intensities = []
    meta = {}
    for line in open(filename):
        line = line.strip()
        if len(line) > 0:
            if ':' in line:
                key, value = line.split(':', 1)
                meta[key] = value.strip()
                if key.strip() == 'Name':
                    name = value.strip()
            else:
                mz, intensity = map(float, line.split(' ', 1))
                mz_values.append(mz)
                intensities.append(intensity)
        else:
            if len(mz_values) >= 6 and len(intensities) > 0:
                index += 1
                spectrum = Spectrum(index, mz_values, intensities)
                spectrum.name = name
                spectrum.meta = meta
                if spectrum.meta['Ion_mode'] == 'P':
                    spectrum.mz_values, spectrum.intensities = remove_peaks(mz_values, intensities)
                    exclude_mz1 = float(meta['PrecursorMZ'])
                    exclude_mz2 = float(meta['ExactMass']) + 1.00784
                    spectrum.mz_values, spectrum.intensities = remove_peaks_P_N(spectrum.mz_values,
                                                                                spectrum.intensities, exclude_mz1)
                    spectrum.mz_values, spectrum.intensities = remove_peaks_P_N(spectrum.mz_values,
                                                                                spectrum.intensities, exclude_mz2)
                    spectra_p.append(spectrum)
                elif spectrum.meta['Ion_mode'] == 'N':
                    spectrum.mz_values, spectrum.intensities = remove_peaks(mz_values, intensities)
                    exclude_mz1 = float(meta['PrecursorMZ'])
                    exclude_mz2 = float(meta['ExactMass']) - 1.00784
                    spectrum.mz_values, spectrum.intensities = remove_peaks_P_N(spectrum.mz_values,
                                                                                spectrum.intensities, exclude_mz1)
                    spectrum.mz_values, spectrum.intensities = remove_peaks_P_N(spectrum.mz_values,
                                                                                spectrum.intensities, exclude_mz2)
                    spectra_n.append(spectrum)

            name = None
            mz_values = []
            intensities = []
            meta = {}

    return spectra_n, spectra_p


def read_file_n_o(mspFile):
    index = 0
    spectra = []
    DB = None
    mz_values = []
    intensities = []
    meta = {}

    for line in open(mspFile):
        line = line.strip()
        if len(line) > 0:
            if ":" in line:
                key, value = line.split(':', 1)
                meta[key] = value
                if key.strip() == "DB#":
                    DB = value.strip().replace("LipidBlast", "1")
            else:
                mz, intensity = map(float, line.split(" ", 1))
                mz_values.append(mz)
                intensities.append(intensity)
        else:
            if DB is not None:
                index += 1
                spectrum = Spectrum(index, mz_values, intensities)
                spectrum.spectrum_id = DB
                spectrum.meta = meta
                spectra.append(spectrum)
            DB = None
            mz_values = []
            intensities = []
            meta = {}
    print("Spectra Collected")
    return spectra


def create_normalized_distance_file(spectra, filename):
    size = len(spectra)
    norms = []
    for spectrum in spectra:
        norms.append(np.sum(spectrum.intensities))

    result_spectra1 = []
    result_spectra2 = []
    result_distance_list = []
    inputs = list()
    for i in range(size):
        inputs.append([i, spectra, norms, size])
    var_pool = mp.Pool()
    print("Start")
    results = var_pool.map(create_normalized_distance_loop, inputs)
    for spectra_1, spectra_2, distance in results:
        result_spectra1 += spectra_1
        result_spectra2 += spectra_2
        result_distance_list += distance
    data_frame_structure = {"SpectrumID 1": result_spectra1, "SpectrumID 2": result_spectra2,
                            "Distance": result_distance_list}
    data_frame_distance = pd.DataFrame(data_frame_structure)
    data_frame_distance.to_csv(filename, index=False)
    n = len(pd.unique(data_frame_distance["SpectrumID 1"]))
    print(n)


def create_normalized_distance_loop(input):
    i, spectra, norms, size = input
    spectra1_list = []
    spectra2_list = []
    distance_list = []
    for j in range(i, size):
        spectrum1 = spectra[i]
        spectrum2 = spectra[j]
        norm1 = norms[i]
        norm2 = norms[j]
        if abs(float(spectrum1.meta["PrecursorMZ"]) - float(spectrum2.meta["PrecursorMZ"])) < 0.01:
            similarity = spectrum1.match(spectrum2, normalized=False)
            similarity /= norm1 * norm2
            distance = 1 - similarity
            if distance <= 0.9:
                spectra1_list.append(spectrum1.spectrum_id)
                spectra2_list.append(spectrum2.spectrum_id)
                distance_list.append(distance)
    return spectra1_list, spectra2_list, distance_list


# create_distance_file(spectra_n, '/Users/ciaraconway/Documents/testing_test_N.csv')
# create_distance_file(spectra_p, '/Users/ciaraconway/Documents/testing_test_P.csv')

# create_normalized_distance_file(spectra_n, '/Users/ciaraconway/Documents/testing_test_Nt.csv')
# create_normalized_distance_file(spectra_p, '/Users/ciaraconway/Documents/testing_test_Pt.csv')

def main():
    # spectra_n, spectra_p = read_file('/Users/ciaraconway/Desktop/Extract.txt')
    # create_normalized_distance_file(spectra_n, '/Users/ciaraconway/Documents/testing_test_Nt.csv')
    # create_normalized_distance_file(spectra_p, '/Users/ciaraconway/Documents/testing_test_Pt.csv')
    msp = "/Users/ciaraconway/Documents/MoNA - Clustering/MoNA-export-LipidBlast-09121-reg_new_331.msp"
    msp1 = "/Users/ciaraconway/Documents/MoNA - Clustering/MoNA-export-LipidBlast-09121-neutral2_neg_new_331.msp"
    spectra = read_file_n_o(msp)
    create_normalized_distance_file(spectra,'/Users/ciaraconway/Documents/MoNA - Clustering/MoNA-export-LipidBlast-09121-reg_new_all_331_mp.csv')
    #spectra1 = read_file_n_o(msp1)
    #create_normalized_distance_file(spectra1,"/Users/ciaraconway/Documents/MoNA - Clustering/MoNA-export-LipidBlast-09121-neutral2_neg_new_all_331_mp.csv")


if __name__ == '__main__':
    main()
