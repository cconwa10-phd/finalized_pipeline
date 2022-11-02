import numpy as np
from sklearn.impute import SimpleImputer
import pandas as pd
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

def read_file(file, type, sheet):
    if type.lower() == "csv":
        df = pd.read_csv(file)
    elif type.lower() == "xlsx":
        df = pd.read_excel(file, sheet_name=sheet)
    return df

def univariate_impute(dataframe, strategy):
    imp = SimpleImputer(missing_values=np.nan, strategy=strategy)
    imp.fit_transform(dataframe)
    return dataframe

def multivariate_impute(dataframe, max_it):
    imp = IterativeImputer(max_iter=max_it, random_state=0)
    imp.fit_transform(dataframe)
    pass
