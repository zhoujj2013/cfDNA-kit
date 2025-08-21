#%%
from sklearn.decomposition import NMF
import numpy as np
from sklearn.linear_model import Lasso
import pandas as pd
import argparse


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


def parse_args():
    parser = argparse.ArgumentParser(description="NMF Matrix Processing")

    parser.add_argument(
        "--expression_matrix",
        type=str,
        required=True,
        help="Path to expression matrix CSV"
    )

    parser.add_argument(
        "--feature_matrix",
        type=str,
        required=True,
        help="Path to feature matrix CSV"
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to weight CSV file"
    )

    return parser.parse_args()

args = parse_args()

expression_matrix_path = args.expression_matrix
feature_matrix_path = args.feature_matrix
csv_file = args.output
expression_matrix_pd = pd.read_csv(expression_matrix_path,index_col=0,header=0)
feature_matrix_pd = pd.read_csv(feature_matrix_path,index_col=0,header=0)
tissues = expression_matrix_pd.index.to_list()
feature = feature_matrix_pd.index.to_list()
expression_matrix = expression_matrix_pd.T.values
feature_matrix = feature_matrix_pd.values
Weight = cal_svd_weight(feature_matrix,expression_matrix)
weight_dict = {}
for i,value in enumerate(Weight):
    weight_dict[tissues[i]] = {}
    for j,v in enumerate(Weight[i]):
        weight_dict[tissues[i]][feature[j]] = v
weight_dict = pd.DataFrame(weight_dict).T
weight_dict.to_csv(csv_file)
#%%