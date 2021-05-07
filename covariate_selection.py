import pandas as pd
from sklearn.linear_model import LogisticRegressionCV, LassoCV
import patsy

def lasso_feature_selection(endpoint, matrix_path):
    matrix = pd.read_csv(matrix_path,
                         dtype={'race_factor': 'category'})
    with pd.option_context('mode.use_inf_as_na', True):
        matrix = matrix.dropna(subset=[endpoint])
    covariates = ["age", "sex", "age_sex", "age_squared", "is_hispanic", "race_factor", "any_comorbidity", "multi_batch"]
    covariates.extend([f"pc{i}" for i in range(1,11)])
    covariates.extend(matrix.columns[matrix.columns.str.startswith("flowcell_")])
    formula_rhs = " + ".join(covariates)
    y, x = patsy.dmatrices(endpoint + "~" + formula_rhs, matrix)

    if pd.api.types.is_bool_dtype(matrix[endpoint]):
        model = LogisticRegressionCV(penalty="l1", solver="liblinear")
    else:
        model = LassoCV()

    model.fit(x, y)

    selected_covars = set()
    for name, coef in zip(x.design_info.column_names, model.coef_.ravel()):
        if coef != 0:
            if name.startswith("pc"):
                # must include all lower PCs
                for i in range(1, int(name[2:]) + 1):
                    selected_covars.add(f"pc{i}")
            elif name.startswith("flowcell_"):
                # include all flowcells if any are selected
                for field in matrix.columns:
                    if field.startswith("flowcell_"):
                        selected_covars.add(field)
            else:
                selected_covars.add(name)
    return {'lasso': list(selected_covars)}
