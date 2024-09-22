import pandas as pd
import numpy as np
np.seterr(divide='ignore', invalid='ignore') # ignores "division by 0" error
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn import metrics
from sklearn.inspection import permutation_importance


def run_ada(bin_mat, pheno_mat, bin_mat_test, pheno_mat_test, max_depth, n_estimators, learning_rate):

    # perpare pheno_mat by converting pheno to string
    pheno_mat = pheno_mat.pheno.values.tolist()
    pheno_mat = [str(x) for x in pheno_mat]

    # prepare bin_mat by dropping converting factor to dichotomous values
    def factor_to_numeric(x):
        if x == 'hit': return -1
        if x == 'miss':   return 1

    bin_mat = bin_mat.applymap(factor_to_numeric)

    model = AdaBoostClassifier(
        DecisionTreeClassifier(max_depth=int(max_depth)),
        n_estimators=int(n_estimators),
        learning_rate=learning_rate,
        algorithm="SAMME.R")

    model.fit(bin_mat, pheno_mat)
    # load the test data
    bin_mat_test = bin_mat_test.applymap(factor_to_numeric)

    pheno_mat_test = pheno_mat_test.pheno.values.tolist()
    pheno_mat_test = [str(x) for x in pheno_mat_test]

    df_probs = pd.DataFrame(model.predict_proba(bin_mat_test))
    df_probs["observed"] = pheno_mat_test

    # analyse the feature importances now
    df_gini = pd.DataFrame({'fature_names': model.feature_names_in_, 'MeanDecreaseGini': model.feature_importances_ })

    return df_probs, df_gini