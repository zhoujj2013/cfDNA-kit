#%%
import os

import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from scipy.fft import fft, fftfreq
from scipy import signal
import numpy as np
import pysam
import argparse
from sklearn.decomposition import NMF
import numpy as np
from sklearn.linear_model import Lasso


def compute_band_mean_powers(wps_arr, fs=1000):
    """
    对 wps_arr 做傅里叶变换，提取 0-250Hz，每 5Hz 一组的均值功率
    :param wps_arr: 一维信号数组
    :param fs: 采样频率，默认 1000Hz（你要根据实际设置）
    :return: 每个 5Hz 频段的功率均值列表
    """
    N = len(wps_arr)
    fft_vals = fft(wps_arr)
    freqs = fftfreq(N, d=1/fs)

    # 取正频率部分
    pos_mask = freqs >= 0
    fft_vals = fft_vals[pos_mask]
    freqs = freqs[pos_mask]

    # 计算功率谱
    power = np.abs(fft_vals) ** 2

    # 提取 0-250Hz 的部分
    freq_mask = freqs <= 250
    freqs = freqs[freq_mask]
    power = power[freq_mask]

    # 分成每 5Hz 一段，计算每段均值
    band_means = []
    for start in range(0, 250, 10):
        band_mask = (freqs >= start) & (freqs < start + 10)
        if np.any(band_mask):
            mean_power = np.mean(power[band_mask])
        else:
            mean_power = 0
        band_means.append(mean_power)

    return band_means
def cal_nmf_weight(input_matrix,tissue_expression_matrix,tissue_count):
    n_components = 10
    # Create the NMF model with specified number of components and random initialization
    nmf = NMF(n_components=n_components, init='random', random_state=0)
    # Fit the NMF model to the input matrix and obtain the factorized matrices
    L = nmf.fit_transform(input_matrix)
    W = nmf.components_
    C = np.zeros((tissue_count, n_components))
    alpha = 0.0001
    for j in range(L.shape[1]):
        y = L[:, j]
        X = tissue_expression_matrix  # Feature matrix
        lasso = Lasso(alpha=alpha, fit_intercept=False, max_iter=200000)
        lasso.fit(X, y)
        C[:, j] = lasso.coef_
    Weight = np.dot(C, W)
    # Reconstruct the original energy matrix using the tissue expression matrix and Weight
    wps_energies_matrix_hat = tissue_expression_matrix @ Weight
    error = np.linalg.norm(input_matrix - wps_energies_matrix_hat) / np.linalg.norm(input_matrix)
    return Weight

#%%

fs = 1000  # 采样率
cutoff_frequency = 100  # 截止频率
order = 10  # 滤波器阶数
b, a = signal.butter(order, cutoff_frequency, fs=fs, btype='low')
# expression_matrix_path = '/mnt/dfc_data2/project/linyusen/project/31_cfdna_wps/project/NMF/tissue_expression_matrix.csv'
# expression_matrix_path = '/mnt/dfc_data2/project/linyusen/project/31_cfdna_wps/project/NMF/tissue_expression_matrix.csv'
# expression_matrix_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/plasma_100/single_cell/single_cell.expression.csv'
expression_matrix_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/nerous/single_cell/single_cell.expression.csv'
expression_matrix_pd = pd.read_csv(expression_matrix_path,index_col=0,header=0)
expression_martrix = expression_matrix_pd.to_dict()
tissues_list = expression_matrix_pd.columns.to_list()
print(expression_matrix_pd.shape)

import pickle

def cal_tissue_matrix(expression_martrix,select_gene):
    tissue_expression_matrix = np.zeros((len(tissues_list), len(uniform_genes)))
    for i, tissue in enumerate(tissues_list):
        for j, feature in enumerate(select_gene):
            tissue_expression_matrix[i, j] = expression_martrix[tissue].get(feature, 0)
    tissue_expression_matrix = tissue_expression_matrix.T
    return tissue_expression_matrix
def cal_feature_matrix(wps_dict,select_gene):
    wps_energies_matrix = []
    for gene in tqdm(select_gene):
        wps_arr = wps_dict[gene]
        wps_arr = signal.filtfilt(b, a, wps_arr)
        band_means = compute_band_mean_powers(wps_arr)  # Compute mean power values in each predefined frequency band from the filtered error signal
        wps_energies_matrix.append(band_means)
    wps_energies_matrix = np.array(wps_energies_matrix)
    return wps_energies_matrix

#%%
import argparse
parser = argparse.ArgumentParser(
        description="Load pickle files, compute results, and save as CSV"
    )
parser.add_argument(
    "--pickle_dir",
    type=str,
    help="输入 pickle 文件或目录路径"
)
parser.add_argument(
    "--csv_file",
    type=str,
    help="结果保存为 CSV 文件路径"
)
parser.add_argument(
    "--scale_score",
    type=float,
    help="相关缩放比例参数（float）"
)

args = parser.parse_args()
pickle_dir = args.pickle_dir
csv_file = args.csv_file
scale_score = args.scale_score



with open(pickle_dir, 'rb') as f:
    gene_wps_dict = pickle.load(f)

for gene in gene_wps_dict:
    gene_wps_dict[gene] = gene_wps_dict[gene]/scale_score
genes = sorted(gene_wps_dict.keys())
uniform_genes = [g for g in genes ]
uniform_genes = list(set(uniform_genes)&set(expression_matrix_pd.index.to_list()))[:5000]

tissue_expression_matrix = cal_tissue_matrix(expression_martrix, uniform_genes)
wps_energies_matrix = cal_feature_matrix(gene_wps_dict, uniform_genes)

Weight = cal_nmf_weight(wps_energies_matrix,tissue_expression_matrix,expression_matrix_pd.shape[1])

weight_dict = {}
tissues = list(expression_martrix.keys())
feature = []
for i in range(25):
    feature.append(i)
for i, value in enumerate(Weight):
    weight_dict[tissues[i]] = {}
    for j, v in enumerate(Weight[i]):
        weight_dict[tissues[i]][feature[j]] = v
weight_dict = pd.DataFrame(weight_dict).T
weight_dict.to_csv(csv_file)

