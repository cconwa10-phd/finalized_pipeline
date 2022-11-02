from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from math import sqrt
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xlrd, openpyxl



def data_cleaning(hfile, hsheet, ifile, isheet):
    ##### ANALYSIS #####
    file = "/Users/cconwa10/Desktop/all_NHANES_KEY.xlsx"
    sheet = "Data_No_Comments"
    df = pd.read_excel(file, sheet_name=sheet)
    hdf = pd.read_excel(hfile, sheet_name=hsheet)
    idf = pd.read_excel(ifile, sheet_name=isheet)
    target = hdf[hdf["target"]].values.to_list
    for i, tar in enumerate(target):

