import numpy as np
from sklearn.impute import SimpleImputer
import pandas as pd
import sklearn
from scipy import stats
import statsmodels.api as sm
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
import seaborn as sns
import matplotlib.pyplot as plt

#Feature_Selection
def univariate_feature_selection(data, target, features):
    best_feat = SelectKBest(score_func=chi2, k=5)
    fit = best_feat.fit(data, target)
    scores = pd.DataFrame(fit.scores_)
    cols = pd.DataFrame(data.columns)
    score_feat = pd.concat([cols, scores], axis=1)
    score_feat.columns = ["Col", "Score"]
    largest = score_feat.nlargest(features, "Score")
    largest = largest.index.values.tolist()
    return largest, score_feat

def corr_heatmap_pearson(dataframe):
    corrmat = dataframe.corr()
    corrmat.to_csv('pearson_corr.csv')
    top_corr_features = corrmat.index
    plt.figure(figsize=(50, 50))
    #g = sns.heatmap(dataframe[top_corr_features].corr(), annot=True, cmap="RdYlGn", fmt='.1g')
    sns.heatmap(dataframe[top_corr_features].corr(), annot=False, cmap="Blues", fmt='.1g', linewidths=1, linecolor='black', robust = True)
    plt.savefig('pearsonR.png')
    #plt.gcf().set_size_inches(15, 8)
    #plt.show()


def corr_heatmap_spearman(dataframe):
    corrmat = dataframe.corr()
    corrmat.to_csv('spearman_corr.csv')
    top_corr_features = corrmat.index
    plt.figure(figsize=(50, 50))
    #g = sns.heatmap(dataframe[top_corr_features].corr(method = 'spearman'), annot=True, cmap="RdYlGn", fmt='.1g')
    sns.heatmap(dataframe[top_corr_features].corr(method = 'spearman'), annot=False, cmap="Blues", fmt='.1g', linewidths=1, linecolor='black', robust = True)
    plt.savefig('spearmanR.png')
    #plt.gcf().set_size_inches(15, 8)
    #plt.show()


###ANALYSIS####

file = "/Users/ciaraconway/Desktop/all_NHANES_KEY.xlsx"
sheet = "Data_No_Comments"
df = pd.read_excel(file, sheet_name=sheet)
df = df[1:len(df)]

corr_heatmap_pearson(df)
corr_heatmap_spearman(df)


