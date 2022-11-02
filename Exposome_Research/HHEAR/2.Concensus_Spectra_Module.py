# Steps for building consensus spectra:
# Collect all m/z peaks from all spectra in a cluster.
# Cluster the m/z values of those peaks with m/z threshold = 0.01, using the complete-linkage hierarchical clustering.
# For each m/z-cluster, add a peak to the consensus spectrum with the following values:
#   m/z = average of all m/z values in the m/z cluster
#   intensity = sum of the intensities of all m/z peaks in the cluster, divided by the number of spectra.
import csv
import pandas as pd
import numpy as np
from spectrum import Spectrum
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
import matplotlib.pyplot as plt
import matplotlib.style as style
import xlsxwriter
style.use('default')


    # return spectra
# for j, labels in enumerate(labels)
def collect_clustered_mz(filename):
    tot_df = pd.DataFrame()
    index = 0
    spectra_collection = []
    data_csv = pd.read_csv(filename)
    labels = data_csv['Cluster'].values
    print(labels)
    labels_list = data_csv['Cluster'].values.tolist()
    clusters = list(dict.fromkeys(labels_list))
    print(clusters)
    for label in clusters:
        index += 1
        mz_cluster = []
        intensity_cluster = []
        cluster_num = []
        precursor_mz = []
        for j in range(0, len(labels)):
            if labels[j] == label:
                spectrum_ids = data_csv.loc[data_csv.Cluster == labels, 'Name'].values[j]
                mz_values = data_csv.loc[data_csv.Cluster == labels, 'mz'].values[j]
                intensities = data_csv.loc[data_csv.Cluster == labels, 'int'].values[j]
                precursor = data_csv.loc[data_csv.Cluster == labels, 'PrecursorMZ'].values[j]
                mz_cluster.append(mz_values)
                intensity_cluster.append(intensities)
                precursor_mz.append(precursor)
                cluster_num.append(spectrum_ids)
        if len(cluster_num) > 1:
            mz_value_str, int_value_str = list_clean_up(mz_cluster, intensity_cluster)
            mz_value_list = list(mz_value_str.split(","))
            int_value_list = list(int_value_str.split(","))
            similarity_matrix_array = normalized_similarity_mz(mz_value_str)
            new_clusters = distance_matrix_clustering(similarity_matrix_array)
            data_frame = {"mz": mz_value_list, "int": int_value_list, "New_Labels": new_clusters}
            df = pd.DataFrame(data_frame, columns=['mz', 'int', 'New_Labels'])
            con_mz, con_int = clustered_mz_int(df, cluster_num)
            stem_graph_cluster(con_mz, con_int, label)
            spectra_new = Spectrum(index, con_mz, con_int)
            spectra_collection.append(spectra_new)
            precursor_mean = np.array(precursor_mz).mean()
            tot = {"Name": cluster_num[0],"mz_c": con_mz, "int_c": con_int, "pre_mean": precursor_mean, "cluster_list": cluster_num}
            tot_df = tot_df.append(tot, ignore_index=True)
        else:
            # mz_cluster = mz_cluster[0]
            # intensity_cluster = intensity_cluster[0]
            # for i, intent in enumerate(intensity_cluster):
            #     intensity_cluster[i] = (float(intent)/max(intensity_cluster))*100
            #stem_graph_cluster(mz_cluster, intensity_cluster, label)
            spectra = Spectrum(index, mz_cluster, intensity_cluster)
            spectra_collection.append(spectra)
    tot_df.to_csv('/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Clusters/UCD_MSMS_QTOF_spectra_data_clustering_conc.csv')
    return spectra_collection, tot_df

def list_clean_up(mz_values,intensities):

    mz_value_str = ','.join(mz_values).replace('[', '').replace(']', '')
    int_value_str = ','.join(intensities).replace('[', '').replace(']', '')
    # cycle_mz = 1
    # cycle_int = 1
    # for mz_value in mz_values:
    #     if cycle_mz == 1:
    #         mz_value = mz_value.strip()
    #         mz_value_str += mz_value
    #         cycle_mz = cycle_mz + 1
    #     else:
    #         mz_value = mz_value.strip()
    #         mz_value_str = mz_value_str + "," + mz_value
    #         cycle_mz = cycle_mz + 1
    # for int_value in intensities:
    #     if cycle_int == 1:
    #         int_value = int_value.strip()
    #         int_str += int_value
    #         cycle_int = cycle_int + 1
    #     else:
    #         int_value = int_value.strip()
    #         int_str = int_str + "," + int_value
    #         cycle_int = cycle_int + 1
    return mz_value_str, int_value_str



def normalized_similarity_mz(mz_values):
    single_mz = mz_values.split(",")
    size = len(single_mz)
    similarity_matrix_mz = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            mz_1 = single_mz[i]
            mz_2 = single_mz[j]
            similarity = abs(np.float64(mz_1) - np.float64(mz_2))
            similarity_matrix_mz[i, j] = similarity
    return similarity_matrix_mz

# def normalized_similarity_mz(mz_values):
#     mz_value = mz_values.split(",")




def distance_matrix_clustering(similarity_matrix_array):
    # cluster = []
    # indices = []
    distance_matrix = similarity_matrix_array
    data_new = ssd.squareform(distance_matrix, force = 'tovector', checks = True)
    Z = sch.linkage(data_new, 'complete')
    #dn = sch.dendrogram(Z)
    hc = sch.fcluster(Z,t = 0.01, criterion = 'distance', depth=2, R=None, monocrit=None)
    cluster = hc.tolist()
    # n_clusters = len(np.unique(hc))
    # n_clusters_loop = np.unique(hc)
    # for n_cluster in n_clusters_loop:
    #     cluster_indices = np.arange(start = 0, stop = len(hc), step = 1)
    #     cluster_group = hc
    #     cluster.append(cluster_group)
    #     indices.append(cluster_indices)
    return cluster


def clustered_mz_int(data_frame, cluster_num):
    con_mz = []
    con_int = []
    dataframe_read = data_frame
    labels = dataframe_read['New_Labels'].values
    # print(labels)
    labels_list = dataframe_read['New_Labels'].values.tolist()
    clusters = list(dict.fromkeys(labels_list))
    # print(clusters)
    for label in clusters:
        new_mz_cluster = []
        new_intensity_cluster = []
        for j in range(0, len(labels)):
            if labels[j] == label:
                mz_values = dataframe_read.loc[dataframe_read.New_Labels == labels, 'mz'].values[j]
                intensities = data_frame.loc[dataframe_read.New_Labels == labels, 'int'].values[j]
                new_mz_cluster.append(float(mz_values))
                new_intensity_cluster.append(float(intensities))
        ave_mz = sum(new_mz_cluster)/len(new_mz_cluster)
        ave_int = sum(new_intensity_cluster)/len(new_intensity_cluster)
        con_mz.append(float(ave_mz))
        con_int.append(float(ave_int))
    # for i, intent in enumerate(con_int):
    #     con_int[i] = (intent/max(con_int))*100
    print(con_mz)
    print(con_int)
    return con_mz,con_int

def stem_graph_cluster(mz, intensity, label):
    x = mz
    y = intensity
    print(x)
    print(y)
    plt.stem(x, y, use_line_collection=True)
    plt.title(label)
    plt.ylabel('intensity')
    plt.xlabel('mass fragmentation')
    #plt.yticks([ 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110])
    plt.show()

def split_value(cell):
    cell = cell.split("[", 1)[1]
    cell = cell.split("]", 1)[0]
    return cell.strip()

def dbID_loop_msp(dbIDs, mspID):
    dbIDString = None
    for i, dbID in enumerate(dbIDs):
        if dbID == mspID:
            dbIDString = mspID
    return dbIDString

def new_msp_file(tot_df, mspFile : str):
    concensus_msp = ""
    concensus_msp_tot = ""
    name = None
    name_ID = tot_df["Name"].values
    for i, ID in enumerate(name_ID):
        for line in open(mspFile):
            line = line.strip()
            if len(line) > 0:
                if ":" in line:
                    key, value = line.split(":", 1)
                    if key.strip() == "Name":
                        if ID.strip() == value.strip():
                            name = value.strip()
                        concensus_msp += line + "\n"
                    elif key.strip() == "PrecursorMZ":
                        pMZ = tot_df.loc[tot_df.Name == name_ID, "pre_mean"].values[i]
                        concensus_msp += "PrecursorMZ: " + str(pMZ) + "\n"
                    else:
                        concensus_msp += line + "\n"
                else:
                    pass
            else:
                if name is not None:
                    for j, k in zip(tot_df.loc[tot_df.Name == name_ID, "mz_c"].values[i], tot_df.loc[tot_df.Name == name_ID, "int_c"].values[i]):
                        concensus_msp += str(j) + "\t" + str(k) + "\n"
                    concensus_msp_tot += concensus_msp + "\n"
                concensus_msp = ""
                name = None
    return concensus_msp_tot


def write_out_msp_new(mspFile, file):
    out_file = open(file, "w")
    out_file.write(mspFile)
    out_file.close()
    print("Outfile: ", file, " has been saved")



def main():
    alignmentID = 14561
    concensus_spectra, tot_df = collect_clustered_mz('/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Clusters/UCD_MSMS_QTOF_spectra_data_clustering_14561.csv')
    mspFile = "/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Files/UCD_MSMS_QTOF_TOT_MW.msp"
    conc_msp_file = new_msp_file(tot_df, mspFile)
    outfile = "/Users/ciaraconway/Documents/MSMS_lipids_allLabs/MSP_Clusters/UCD_MSMS_QTOF_14561_concensus.msp"
    write_out_msp_new(conc_msp_file, outfile)
if __name__ == '__main__':
    main()
