from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
import pandas as pd
import numpy as np
import math

def log_reg(bin_mat, pheno_mat, bin_mat_test, pheno_mat_test):
    # bin_mat = pd.read_csv("C:/Users/bujdo/OneDrive/Dokumenty/PhD_OToole/aurora/bin_mat.csv")
    # pheno_mat = pd.read_csv("C:/Users/bujdo/OneDrive/Dokumenty/PhD_OToole/aurora/pheno_mat.csv")
    # bin_mat_test = pd.read_csv("C:/Users/bujdo/OneDrive/Dokumenty/PhD_OToole/aurora/bin_mat_test.csv")
    # pheno_mat_test = pd.read_csv("C:/Users/bujdo/OneDrive/Dokumenty/PhD_OToole/aurora/pheno_mat_test.csv")

    pheno_mat = [str(x) for x in pheno_mat["pheno"]]
    pheno_mat_test = [str(x) for x in pheno_mat_test["pheno"]]

    def factor_to_numeric(x):
        if x == 'hit': return -1
        if x == 'miss': return 1
        else: return x



    model = LogisticRegression(max_iter = 1000, penalty = "l1", C = 0.1, solver = "liblinear")
    # fit model
    model.fit(bin_mat, pheno_mat)

    # confusion_matrix(model.predict(bin_mat_test), pheno_mat_test)
    # save results
    df_probs = pd.DataFrame(model.predict_proba(bin_mat_test))
    df_probs["observed"] = pheno_mat_test
    return df_probs



