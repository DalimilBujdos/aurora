import pandas as pd
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn import metrics


def run_ada_miss(bin_mat_ada, pheno_mat_ada, bin_mat_test_ada, pheno_mat_test_ada, max_depth, n_estimators, learning_rate):
    # perpare pheno_mat by converting pheno to string
    pheno_mat = pheno_mat_ada.pheno.values.tolist()
    pheno_mat = [str(x) for x in pheno_mat]

    # prepare bin_mat by dropping converting factor to dichotomous values
    def factor_to_numeric(x):
        if x == 'hit': return -1
        if x == 'miss':   return 1

    bin_mat = bin_mat_ada.applymap(factor_to_numeric)

    model = AdaBoostClassifier(
        DecisionTreeClassifier(max_depth=int(max_depth)),
        n_estimators=int(n_estimators),
        learning_rate=learning_rate,
        algorithm="SAMME.R")

    model.fit(bin_mat, pheno_mat)
    # load the test data
    bin_mat_test = bin_mat_test_ada.applymap(factor_to_numeric)

    df_probs = pd.DataFrame(model.predict_proba(bin_mat_test))
    return df_probs
