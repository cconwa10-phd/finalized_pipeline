import numpy as np
from sklearn.impute import SimpleImputer
import pandas as pd
import sklearn
import sklearn.pipeline as sp
from scipy import stats
import statsmodels.api as sm
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
import seaborn as sns
import matplotlib.pyplot as plt
import xlrd, openpyxl

#Optional Imputation
def univariate_impute(dataframe, strategy):
    imp = SimpleImputer(missing_values=np.nan, strategy=strategy)
    return imp.fit_transform(dataframe)


#Looped Linear Regression
def lin_reg(X , y):
    reg = sklearn.linear_model.LinearRegression().fit(X, y)
    return [reg.score(X, y), reg.coef_, reg.intercept_, reg.predict(X)]

def p_value_LR_scipy(X, y):
    lm = sklearn.linear_model.LinearRegression()
    lm.fit(X, y)
    params = np.append(lm.intercept_, lm.coef_)
    pred = lm.predict(X)

    p_X = X.to_numpy()
    p_y = y.to_numpy()

    n_X = np.append(np.ones((len(X), 1)), X, axis=1)
    mse = (sum((p_y-pred)**2))/(len(n_X) - len(n_X[0]))
    try:
        vb = mse * (np.linalg.inv(np.dot(n_X.T, n_X)).diagonal())
        sb = np.sqrt(vb)
        tb = params / sb
        pval = [2 * (1 - stats.t.cdf(np.abs(i), (len(n_X) - len(n_X[0])))) for i in tb]
        return pval
    except:
        return 0
    # vb = mse*(np.linalg.inv(np.dot(n_X.T, n_X)).diagonal())
    # sb = np.sqrt(vb)
    # tb = params/sb
    #
    # pval = [2*(1-stats.t.cdf(np.abs(i),(len(n_X)-len(n_X[0])))) for i in tb]
    # #pval = np.round(pval, 10)
    # return pval

def p_val_LR_stats(X,y):
    lm = sm.OLS(y,X)
    lmf = lm.fit()
    pval = lmf.summary2().tables[1]['P>|t|']
    return pval


def lin_reg_loop_iter(df):
    outer_n = []
    inner_n = []
    lm_score = []
    pval_sci = []
    pval_stat = []
    for (coln, cold) in df.iteritems():
        for (coln1, cold1) in df.iteritems():
            outer_n.append(coln)
            inner_n.append(coln1)
            lin = lin_reg(cold.values, cold1.values)
            pval1 = p_value_LR_scipy(cold.values, cold1.values)
            pval2 = p_val_LR_stats(cold.values, cold1.values)
            lm_score.append(lin[0])
            pval_sci.append(pval1)
            pval_stat.append(pval2)
    data = list(zip(outer_n,inner_n,lm_score,pval_sci,pval_stat))
    final_df = pd.DataFrame(data, columns=["ID1", "ID2", "LinearModScore", "Pval_Scipy", "Pval_Stat"])
    return final_df

def lin_reg_loop(df):
    feat = df.columns.values.tolist()
    feat = feat[1:len(feat)]
    outer_n = []
    inner_n = []
    lm_score = []
    pval_sci = []
    pval_stat = []
    for i, cold1 in enumerate(feat):
        for j, cold2 in enumerate(feat):
            n_df = df[[cold1, cold2]].copy()
            n_df = n_df.dropna()
            X = n_df[[cold1]]
            y = n_df[[cold2]]
            outer_n.append(cold1)
            inner_n.append(cold2)
            lin = lin_reg(X, y)
            pval1 = p_value_LR_scipy(X, y)
            pval2 = p_val_LR_stats(X, y)
            lm_score.append(lin[0])
            pval_sci.append(pval1)
            pval_stat.append(pval2)
            print(str(lin[0])+" " +str(pval1) + " "+str(pval2))
    data = list(zip(outer_n,inner_n,lm_score,pval_sci,pval_stat))
    final_df = pd.DataFrame(data, columns=["ID1", "ID2", "LinearModScore", "Pval_Scipy", "Pval_Stat"])
    return final_df


###Other Regressions###
def lin_reg_other(X , y):
    reg = sklearn.linear_model.LinearRegression().fit(X, y)
    return reg.score(X,y)
def poly_reg(X, y, degree, model):
    poly = sklearn.preprocessing.PolynomialFeatures(degree=degree)
    X_poly = poly.fit_transform(X)
    if model.lower() == "log":
        model = sklearn.linear_model.LogisticRegression()
        model.fit(X_poly, y)
        score = model.score(X_poly, y)
        return score
    elif model.lower() == "lin":
        model = sklearn.linear_model.LinearRegression()
        model.fit(X_poly, y)
        score = model.score(X_poly, y)
        return score
def rid_reg(X, y):
    reg = sklearn.linear_model.Ridge(alpha=1.0)
    reg.fit(X, y)
    return reg.score(X,y)
def las_reg(X, y):
    reg = sklearn.linear_model.Lasso(alpha=0.1)
    reg.fit(X, y)
    path = reg.path(X, y)
    score = reg.score(X, y)
    return score
def ela_reg(X,y):
    reg = sklearn.linear_model.ElasticNet(random_state=0)
    reg.fit(X, y)
    path = reg.path(X, y)
    score = reg.score(X, y)
    return score
def sv_reg(X, y):
    reg = sp.make_pipeline(sklearn.preprocessing.StandardScaler(), sklearn.svm.SVR(C = 1.0, epsilon=0.2))
    reg.fit(X,y)
    return reg.score(X,y)
def bay_reg(X, y):
    reg = sklearn.linear_model.BayesianRidge()
    reg.fit(X, y)
    return reg.score(X,y)
def lar_reg(X, y):
    reg = sklearn.linear_model.LassoLars(alpha=0.1,normalize=False)
    reg.fit(X, y)
    return reg.score(X,y)
def reg_loop(df):
    feat = df.columns.values.tolist()
    feat = feat[1:len(feat)]
    select_feat = ['LBXNFOA','LBXNFOS','LBXPFDE','LBXPFHS','LBXPFNA','LBXPFUA','LBXBFOA', 'LBXMFOS', 'LBXMPAH']
    select_feat_1 = ['LBXBFOA', 'LBXMFOS', 'LBXMPAH']
    data = []
    for i, cold1 in enumerate(feat):
        for j, cold2 in enumerate(feat):
            n_df = pd.DataFrame(df[[cold1,cold2]].copy())
            n_df = n_df.dropna()
            X = n_df[[cold1]]
            y = n_df[[cold2]].values.ravel()
            #X = n_df.iloc[:, 0].to_numpy().reshape(1, -1)
            #y = n_df.iloc[:, 1].to_numpy().reshape(1, -1)
            if len(X) != 0 and len(y) != 0 and str(cold1) != str(cold2):
                lin_s = lin_reg(X,y)
                pol1_lin_s = poly_reg(X, y, 1, "lin")
                pol2_lin_s = poly_reg(X, y, 2, "lin")
                pol3_lin_s = poly_reg(X, y, 3, "lin")
                rid_s = rid_reg(X,y)
                las_s = las_reg(X,y)
                ela_s = ela_reg(X,y)
                sv_s = sv_reg(X,y)
                bay_s = bay_reg(X,y)
                lar_s = lar_reg(X,y)
                data.append((cold1,cold2,lin_s,pol1_lin_s,pol2_lin_s,pol3_lin_s,rid_s,las_s,ela_s,sv_s,bay_s,lar_s))
                print(cold1,cold2,lin_s,pol1_lin_s,pol2_lin_s,pol3_lin_s,rid_s,las_s,ela_s,sv_s,bay_s,lar_s)
            else:
                print("Array is 0")
    final_df = pd.DataFrame(data, columns=['cold1','cold2','lin_s','pol1_lin_s','pol2_lin_s','pol3_lin_s','rid_s','las_s','ela_s',
                                           'sv_s', 'bay_s', 'lar_s'])
    return final_df

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

def corr_heatmap(dataframe):
    corrmat = dataframe.corr()
    top_corr_features = corrmat.index
    plt.figure(figsize=(20, 20))
    #g = sns.heatmap(dataframe[top_corr_features].corr(), annot=True, cmap="RdYlGn")
    sns.heatmap(dataframe[top_corr_features].corr(), annot=False, cmap="RdYlGn")
    plt.show()


#Looping Feature Selection Target Column
def feature_select_loop(dataframe, key):
    top_features = []
    target_list = []
    for i in key:
        data = dataframe.drop(["SEQN"], 1)
        data = data.drop(i, 1)
        target = dataframe[i].values
        largest, score_feat = univariate_feature_selection(data, target, 10)
        top_features.append(largest)
        target_list.append(i)
    return top_features, target_list


##### ANALYSIS #####
file = "/Users/cconwa10/Desktop/all_NHANES_KEY.xlsx"
sheet = "Data_No_Comments"
df = pd.read_excel(file, sheet_name=sheet)
df.set_index('SEQN')

impute_data = univariate_impute(df, "mean")
df_new = pd.DataFrame(impute_data, columns=df.columns.values.tolist())
df_new.set_index('SEQN')

# corr_heatmap(df)
# corr_heatmap(df_new)

#lin_df = lin_reg_loop(df)
#lin_df.to_csv("/Users/ciaraconway/Documents/all_databases/NHANES_data/2017-2018/laboratory_data/all_NHANES_KEY_linreg_nonimp.csv")

reg_df = reg_loop(df)
reg_df.to_csv("/Users/cconwa10/Desktop/all_NHANES_KEY_reg_scores_nonimp_all.csv")
