#%%
from sklearn.decomposition import NMF
import numpy as np
from sklearn.linear_model import Lasso

def cal_svd_weight(input_matrix,tissue_expression_matrix):
    U, S, VT = np.linalg.svd(input_matrix, full_matrices=False)
    Sigma = np.diag(S)
    US = U @ Sigma

    L = US
    W = VT
    C = np.zeros((tissue_expression_matrix.shape[1], W.shape[1]))
    alpha = 0.0001  # Lasso 正则项系数，可以调节或使用 LassoCV 自动选取
    for j in range(L.shape[1]):
        y = L[:, j]  # 第 j 列目标向量 (857,)
        X = tissue_expression_matrix  # 特征矩阵 (857, 24)
        lasso = Lasso(alpha=alpha, fit_intercept=False, max_iter=200000)
        lasso.fit(X, y)
        C[:, j] = lasso.coef_  # 存入第 j 列
    Weight = np.dot(C, W)
    # 拟合结果 B @ C 应接近 A
    wps_energies_matrix_hat = tissue_expression_matrix @ Weight
    error = np.linalg.norm(input_matrix - wps_energies_matrix_hat) / np.linalg.norm(input_matrix)
    return Weight
