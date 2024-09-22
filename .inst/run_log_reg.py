import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
np.seterr(divide='ignore', invalid='ignore') # ignores "division by 0" error

def run_log_reg(bin_mat_log_reg, pheno_mat_log_reg, bin_mat_test_log_reg, pheno_mat_test_log_reg, C_val):
    # perpare pheno_mat by converting pheno to string
    pheno_mat = pheno_mat_log_reg.pheno.values.tolist()
    pheno_mat = [str(x) for x in pheno_mat]

    # prepare bin_mat by dropping converting factor to dichotomous values
    def factor_to_numeric(x):
        if x == 'hit': return -1
        if x == 'miss':   return 1

    bin_mat = bin_mat_log_reg.applymap(factor_to_numeric)

    model = LogisticRegression(max_iter=10000,
                               penalty="l1",
                               C=C_val,
                               solver="liblinear",
                               multi_class="ovr")

    model.fit(bin_mat, pheno_mat)
    # load the test data
    bin_mat_test = bin_mat_test_log_reg.applymap(factor_to_numeric)

    pheno_mat_test = pheno_mat_test_log_reg.pheno.values.tolist()
    pheno_mat_test = [str(x) for x in pheno_mat_test]

    df_probs = pd.DataFrame(model.predict_proba(bin_mat_test))
    df_probs["observed"] = pheno_mat_test

    # get model coefficients
    coefs = np.transpose(model.coef_)
    col_nms = [str(x) for x in range(1,coefs.shape[1]+1)]
    df_coefs = pd.DataFrame(coefs, columns = col_nms)

    return df_probs, df_coefs