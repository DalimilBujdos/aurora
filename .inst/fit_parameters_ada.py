import pandas as pd
import numpy as np
import statistics
from random import sample
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn import metrics

def fit_parameters_ada(bags, predictors, pheno_mat, phenotypes, repeats):
    phenotypes = [str(x) for x in phenotypes]
    repeats = int(repeats)

    pheno_mat["pheno"] = [str(x) for x in pheno_mat["pheno"]]


    def factor_to_numeric(x):
        if x == 'hit': return -1
        if x == 'miss':   return 1
        else : return x

    predictors = predictors.applymap(factor_to_numeric)

    max_depth_arg = 1 # incresing this only decreases the prob of missclassified in their original class
    n_estimators_arg = [10,50,100,500]
    n_estimators_arg = [i for i in n_estimators_arg if i < len(predictors.axes[1])]
    learning_rate_arg = [0.01, 0.1, 0.5, 1, 1.5]
    if type(n_estimators_arg) == list:
        len_n_estimators_arg = len(n_estimators_arg)
    else:
        len_n_estimators_arg = 1

    nrows = max_depth_arg * len_n_estimators_arg * len(learning_rate_arg)
    Grid = pd.DataFrame(np.array(np.meshgrid(max_depth_arg, n_estimators_arg, learning_rate_arg)).reshape(3,nrows).T)

    def run_ada_param(x, train_data, train_data_pheno, test_data, test_data_pheno, phenotypes):
        model = AdaBoostClassifier(
                DecisionTreeClassifier(max_depth=int(x[0])),
                n_estimators = int(x[1]),
                learning_rate= x[2],
                algorithm="SAMME.R")

        model.fit(train_data, train_data_pheno)
        predictions = pd.DataFrame(model.predict_proba(test_data))
        def get_pred_pheno(y, phenotypes):
            y = np.array(y.tolist())
            return phenotypes[y.argmax()]

        df_tmp = predictions.iloc[:, range(len(phenotypes))]
        predictions["predicted"] = df_tmp.apply(get_pred_pheno, axis=1, args=[phenotypes])
        predictions["observed"] = test_data_pheno
        def  accu_calc(y, phenotypes):
            y = y.tolist()
            if y[len(phenotypes)] == y[len(phenotypes) + 1]: return 1
            else : return 0

        correct_predictions = np.sum(predictions.apply(accu_calc, axis=1, args=[phenotypes]).to_numpy())
        accuracy = correct_predictions/len(predictions)
        return accuracy


    df_final = pd.DataFrame()
    for i in range(repeats):
        bags_list = bags.iloc[:, i].tolist()
        pheno_mat_ids = pheno_mat.iloc[:, 0].tolist()
        pheno_mat_pheno = pheno_mat.iloc[:, 1].tolist()

        index_select = []
        for j in range(len(bags_list)):
            index = pheno_mat_ids.index(bags_list[j])
            index_select.append(index)

        pheno_mat_model = pd.DataFrame(
            dict(ids = pheno_mat.iloc[index_select, 0], pheno = pheno_mat.iloc[index_select, 1])).reset_index()
        del pheno_mat_model["index"]

        predictors_list = predictors.loc[:, "ids"].tolist()
        index_select = []
        for j in range(len(bags_list)):
            index = predictors_list.index(bags_list[j])
            index_select.append(index)

        predictors_train = predictors.iloc[index_select, :]
        del predictors_train["ids"]
        predictors_train = predictors_train.assign(pheno = pheno_mat_model.iloc[:, 1].tolist())

        predictors_learn = pd.DataFrame()
        predictors_test = pd.DataFrame()
        for j in range(len(phenotypes)):
            subset = predictors_train[predictors_train["pheno"] == phenotypes[j]]
            subset.reset_index(inplace = True)
            del subset["index"]
            index = sample(range(len(subset)), round(0.3*len(subset)))
            predictors_learn = pd.concat([predictors_learn, subset.iloc[index,:]])
            subset_test = subset.drop(index = index, axis = 0)
            predictors_test = pd.concat([predictors_test, subset_test])

        train_data = predictors_learn.drop("pheno", axis=1)
        train_data_pheno = [str(x) for x in predictors_learn["pheno"].tolist()]
        test_data = predictors_test.drop("pheno", axis=1)
        test_data_pheno = [str(x) for x in predictors_test["pheno"].tolist()]

        accuracies = Grid.apply(run_ada_param, axis=1,
                                train_data = train_data,
                                train_data_pheno = train_data_pheno,
                                test_data = test_data,
                                test_data_pheno = test_data_pheno,
                                phenotypes = phenotypes)

        df_final = pd.concat([df_final, accuracies], axis=1)

    df_final["median"] = df_final.apply(statistics.median, axis = 1)
    return_vals = Grid.iloc[df_final["median"].to_numpy().argmax(), :].to_numpy()

    return return_vals

