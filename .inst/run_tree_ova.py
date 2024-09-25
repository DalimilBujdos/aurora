import pandas as pd
import numpy as np
import statistics
from random import sample
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn import metrics

def run_tree_ova(bin_mat_ada, pheno_mat_ada, bin_mat_test_ada, pheno_mat_test_ada, max_depth):
    # perpare pheno_mat by converting pheno to string
    pheno_mat = pheno_mat_ada.pheno.values.tolist()
    pheno_mat = [str(x) for x in pheno_mat]

    # prepare bin_mat by dropping converting factor to dichotomous values
    def factor_to_numeric(x):
        if x == 'hit': return -1
        if x == 'miss':   return 1

    bin_mat = bin_mat_ada.applymap(factor_to_numeric)

    model = DecisionTreeClassifier(max_depth=int(max_depth))
    # define ovo strategy
    ovo_tree = OneVsRestClassifier(model)
    ovo_tree.fit(bin_mat, pheno_mat)
    # load the test data
    bin_mat_test = bin_mat_test_ada.applymap(factor_to_numeric)

    pheno_mat_test = pheno_mat_test_ada.pheno.values.tolist()
    pheno_mat_test = [str(x) for x in pheno_mat_test]

    df_probs = pd.DataFrame(ovo_tree.predict_proba(bin_mat_test))
    return df_probs