import pandas as pd
import math
import time
from sklearn.metrics.cluster import fowlkes_mallows_score

def PPV(pos, fal_pos):
    if pos == 0:
        return 0
    else:
        return pos/(pos + fal_pos)

def FDR(pos, fal_pos):
    if fal_pos == 0:
        return 0
    else:
        return fal_pos/(fal_pos + pos)

def accuracy(pos, neg):
    if pos == 0:
        return 0
    else:
        return pos/(neg+pos)



def info_read_res_tar(target, results):
    align_id = []
    softw = []
    fdr = []
    accu = []
    ppv = []
    pos_mean_score = []
    neg_mean_score = []
    overall_mean_score = []
    tardf = pd.read_csv(target)
    resdf = pd.read_csv(results)
    ver = tardf["Version"].to_list()
    clu = tardf["clusters"].tolist()
    tarpairs = dict(zip(ver, clu))
    count = 0
    for key, value in tarpairs.items():
        software = resdf["Type"].unique().tolist()
        df_tar = resdf[(resdf.Version == key) & (resdf.clusters == value)]
        for soft in software:
            df_tar_ = df_tar[resdf.Type == soft]
            df_tar_m = df_tar_["Score"].mean()
            df_tar_ind = len(df_tar_.index)
            df_result = resdf[(resdf.Type == soft) & (resdf.Version == key) & (resdf.clusters != value)]
            df_result_m = df_result["Score"].mean()
            df_result_all = resdf[(resdf.Type == soft) & (resdf.Version == key)]
            df_result_all_m = df_result_all["Score"].mean()
            df_res_ind = len(df_result.index)
            df_res_all_ind = len(df_result_all.index)
            pos_mean_score.append(df_tar_m)
            neg_mean_score.append(df_result_m)
            overall_mean_score.append(df_result_all_m)
            fdr.append(FDR(df_tar_ind, df_res_ind))
            accu.append(accuracy(df_tar_ind, df_res_ind))
            align_id.append(key)
            softw.append(soft)
            count+=1
            print(count)
    df = pd.DataFrame(data=zip(align_id,softw,fdr,accu,pos_mean_score,neg_mean_score,overall_mean_score), columns=[
        "version", "type", "fdr", "accuracy", "mean score (pos)", "mean score (neg)", "mean score (overall)"
    ])
    return df

def fm_equation(pos, neg):
    if pos == 0:
        return 0
    else:
        return math.sqrt(pos/(pos+neg))

def fm_analysis(target, results):
    df_e = pd.DataFrame()
    tardf = pd.read_csv(target)
    resdf1 = pd.read_csv(results)
    ver = tardf["Version"].to_list()
    clu = tardf["clusters"].tolist()
    tols = [0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
    tarpairs = dict(zip(ver, clu))
    count = 0
    print(len(resdf1))
    for i in range(0, len(resdf1.columns), 3):
        print(i)
        resdf = resdf1.iloc[:,i:i+3]
        cols = list(resdf.columns)
        software = resdf.iloc[:,1].unique().tolist()
        softw = []
        ver = []
        fm_score = []
        for key, value in tarpairs.items():
            df_tar = resdf[(resdf[cols[0]] == key) & (resdf[cols[2]] == value)]
            for soft in software:
                df_tar_ = df_tar[df_tar[cols[1]] == soft]
                df_tar_ind = len(df_tar_.index)
                df_result = resdf[(resdf[cols[1]] == soft) & (resdf[cols[0]] == key) & (resdf[cols[2]] != value)]
                df_res_ind = len(df_result.index)
                #true_cluster = resdf1[(resdf1.iloc[:,i+1] == soft) & (resdf1.iloc[:,i+0] == key)]
                #true_cluster = true_cluster.iloc[:,32]
                exper_cluster = resdf[(resdf[cols[1]] == soft) & (resdf[cols[0]] == key)]
                exper_cluster = exper_cluster[cols[2]].to_list()
                true_cluster = [value]*len(exper_cluster)
                fm = fowlkes_mallows_score(true_cluster, exper_cluster)
                print(fm)
                #fm = fm_equation(df_tar_ind, df_res_ind)
                softw.append(soft)
                ver.append(key)
                fm_score.append(fm)
        df_e[str(tols[count]) + "v"] = ver
        df_e[str(tols[count]) + "s"] = softw
        df_e[str(tols[count]) + "fm"] = fm_score
        count += 1
    return df_e




target = "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/HHEAR_Results/gather_target_matches.csv"
results = "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/HHEAR_Results/all_results_clusters.csv"
#df = info_read_res_tar(target, results)
#df.to_csv("/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/HHEAR_Results/stats_HHEAR.csv")

results2 = "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/HHEAR_Results/fm_hhear_data.csv"

fmr = fm_analysis(target, results2)
fmr.to_csv("/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/HHEAR_Results/fm_hhear_results3.csv")
