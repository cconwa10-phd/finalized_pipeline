import sklearn
import numpy as np


# 1. Linear Regression
def lin_reg(X , y, X_p, predict=True):
    reg = sklearn.linear_model.LinearRegression().fit(X, y)
    if predict:
        return [reg.score(X, y), reg.coef_, reg.intercept_, reg.predict(X_p)]
    else:
        return [reg.score(X, y), reg.coef_, reg.intercept_]
# 2. Polynomial Regression
def poly_reg(X, y, X_p, y_p, degree, model, predict=True):
    poly = sklearn.preprocessing.PolynomialFeatures(degree=degree)
    X_poly = poly.fit_transform(X)
    if model.lower() == "log":
        model = sklearn.linear_model.LogisticRegression()
        model.fit(X_poly, y)
        if predict:
            y_pred = sklearn.linear_model.LogisticRegression().predict(X_p)
            rmse = np.sqrt(sklearn.metrics.mean_squared_error(y_p, y_pred))
            r2 = sklearn.metrics.r2_score(y_p, y_pred)
            return [y_pred, rmse, r2]
    elif model.lower() == "lin":
        model = sklearn.linear_model.LinearRegression()
        model.fit(X_poly, y)
        score = model.score(X_poly, y)
        if predict:
            y_pred = sklearn.linear_model.LinearRegression().predict(X_p)
            rmse = np.sqrt(sklearn.metrics.mean_squared_error(y_p, y_pred))
            r2 = sklearn.metrics.r2_score(y_p, y_pred)
            return [y_pred, rmse, r2, score]

# 3. Logistic Regression
def log_reg(X, y, X_p, y_p):
    reg = sklearn.linear_model.LogisticRegression(random_state=0).fit(X,y)
    reg.fit(X, y)
    predictions = reg.predict(X_p)
    score = reg.score(X_p, y_p)
    cm = sklearn.metrics.confusion_matrix(y_p, predictions)
    return [y_p, predictions, score, cm]

# 4. Quantile Regression
def quart_reg(X, y, quartiles):
    predictions = {}
    y_mean = np.mean(y)
    out_bounds = np.zeros_like(y_mean, dtype=np.bool)
    for quartile in quartiles:
        reg = sklearn.linear_model.QuantileRegressor(quartile=quartile, alpha=0)
        y_p = reg.fit(X, y).predict(X)
        predictions[quartile] = y_p
        if quartile == min(quartiles):
            out_bounds = np.logical_or(out_bounds, y_p >= y)
        elif quartile == max(quartiles):
            out_bounds = np.logical_or(out_bounds, y_p <= y)
    return predictions

# 5. Ridge Regression
def rid_reg(X, y, X_p, y_p):
    reg = sklearn.linear_model.Ridge(alpha=1.0)
    reg.fit(X, y)
    return [reg.score(X,y), reg.predict(X_p)]

# 6. Lasso Regression
def las_reg(X, y, X_p, y_p):
    reg = sklearn.linear_model.Lasso(alpha=0.1)
    reg.fit(X, y)
    path = reg.path(X, y)
    score = reg.score(X, y)
    y_pred = reg.predict(X_p)
    return [path, score, y_pred, y_p]
# 7. Elastic Net Regression
def ela_reg(X,y, X_p, y_p):
    reg = sklearn.linear_model.ElasticNet(random_state=0)
    reg.fit(X, y)
    path = reg.path(X, y)
    score = reg.score(X, y)
    y_pred = reg.predict(X_p)
    return [path, score, y_pred, y_p]
# 8. Principal Components Regression (PCR)
def pc_reg(X, y, X_p, y_p, components):
    pcr = sklearn.pipeline.make_pipeline(sklearn.preprocessing.StandardScaler(), sklearn.decomposition.PCA(n_components=components), sklearn.linear_model.LinearRegression())
    pcr.fit(X,y)
    pca = pcr.named_steps["pca"]
    return [pca.transform(X_p), pcr.predict(X_p), y_p, pca.components_, pca.explained_variance, pcr.score(X_p, y_p)]

# 9. Partial Least Squares (PLS) Regression
def pls_reg(X, y, X_p, y_p, components):
    pls = sklearn.cross_decomposition.PLSRegression(n_components=components)
    pls.fit(X, y)
    return [pls.transform(X_p), pls.predict(X_p), pls.score(X_p, y_p)]
# 10. Support Vector Regression
def sv_reg(X, y, X_p, y_p):
    reg = sklearn.pipeline.make_pipeline(sklearn.preprocessing.StandardScaler(), sklearn.svm.SVR(C = 1.0, epsilon=0.2))
    reg.fit(X,y)
    return [reg.score(X, y), reg.predict(X_p), y_p]
# 11. Ordinal Regression
# 12. Poisson Regression
def pos_reg(X, y, X_p, y_p):
    reg = sklearn.linear_model.PoissonRegressor()
    reg.fit(X, y)
    return [reg.score(X, y), reg.coef_, reg.intercept_, reg.predict(X_p), y_p]
# 13. Negative Binomial Regression
# 14. Quasi Poisson Regression
# 15. Cox Regression
# 16. Tobit Regression
