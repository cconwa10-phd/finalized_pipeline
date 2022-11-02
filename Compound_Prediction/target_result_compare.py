import pandas as pd


def compare_target(target, results):
    emptydf = pd.DataFrame()
    tardf = pd.read_csv(target)
    resdf = pd.read_csv(results)
    ver = tardf["Version"].to_list()
    clu = tardf["clusters"].tolist()
    tarpairs = dict(zip(ver,clu))
    for key, value in tarpairs.items():
        df = resdf[(resdf.Version == key) & (resdf.clusters == value)]
        emptydf = emptydf.append(df)
    return emptydf


target = "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/all_results/gather_target_matches.csv"
results = "/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/all_results/all_results_clusters.csv"
compare_target(target,results).to_csv("/Users/ciaraconway/Documents/all_databases/Spectra/LCMSMS/all_results/tar_clusters_only.csv")

