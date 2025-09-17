#%%
import os
import csv
import pandas as pd
from sklearn.decomposition import NMF
import numpy as np
from sklearn.linear_model import Lasso
import argparse


parser = argparse.ArgumentParser(description="Process gene features and generate result.")
parser.add_argument("--expression_matrix", required=True, help="Path to expression matrix CSV file")
parser.add_argument("--feature_list", required=True, help="Path to feature list file")
parser.add_argument("--csv_file", required=True, help="Path to output result file")
args = parser.parse_args()

expression_matrix_path = args.expression_matrix
feature_file = args.feature_list
csv_file = args.csv_file
# expression_matrix_path = '/mnt/dfc_data1/home/linyusen/31_cfdna_wps/project/NMF/tissue_expression_matrix_specific.csv'
# feature_file = '/mnt/dfc_data2/project/linyusen/project/31_cfdna_wps/project/NMF/multi_feature/feature.lst'
# csv_file = '/mnt/dfc_data2/project/linyusen/project/31_cfdna_wps/project/NMF/multi_feature/tissue_weight.result'

def cal_tissue_matrix(expression_martrix,tissues_list,select_gene):
    tissue_expression_matrix = np.zeros((len(tissues_list), len(select_gene)))
    for i, tissue in enumerate(tissues_list):
        for j, feature in enumerate(select_gene):
            tissue_expression_matrix[i, j] = expression_martrix[tissue].get(feature, 0)
    tissue_expression_matrix = tissue_expression_matrix.T
    return tissue_expression_matrix
def cal_nmf_weight(input_matrix,tissue_expression_matrix):
    n_components = 10
    # Create the NMF model with specified number of components and random initialization
    nmf = NMF(n_components=n_components, init='random', random_state=0)
    # Fit the NMF model to the input matrix and obtain the factorized matrices
    L = nmf.fit_transform(input_matrix)
    W = nmf.components_
    C = np.zeros((24, n_components))
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

expression_matrix_pd = pd.read_csv(expression_matrix_path,index_col=0,header=0)
expression_martrix = expression_matrix_pd.to_dict()
tissues_list = expression_matrix_pd.columns.to_list()
uniform_gene = []
feature_dict_temp = {}
feature_select = {}
f = open(feature_file,'r')
for line in f.readlines():
    select_feature = []
    line = line.strip().split(' ')
    try:
        feature_name = line[0]
        file_path = line[1]
        select_feature = line[2]
        select_feature = select_feature.split(',')
        select_feature = [int(j) for j in select_feature]
    except:
        feature_name = line[0]
        file_path = line[1]
        select_feature = [0]
    feature_dict_temp[feature_name] = {}
    gene_temp = []
    with open(file_path) as f:
        r = csv.reader(f, delimiter='\t')
        for i in r:
            gene = i[0]
            gene_temp.append(gene)
            score = i[1:]
            score = [float(j) for j in score]
            feature_dict_temp[feature_name][gene] = score
            feature_select[feature_name] = select_feature
    uniform_gene.append(gene_temp)
#%%
feature_dict = {}
for feature_name in feature_dict_temp.keys():
    for gene in feature_dict_temp[feature_name]:
        if len(feature_dict_temp[feature_name][gene]) == 1:
            new_name = feature_name
            if new_name not in feature_dict:
                feature_dict[new_name] = {}
            feature_dict[new_name][gene] = feature_dict_temp[feature_name][gene][0]
        else:
            for index,value in enumerate(feature_dict_temp[feature_name][gene]):
                if index not in feature_select[feature_name]:
                    continue
                new_name = feature_name + '_' + str(index)
                if new_name not in feature_dict:
                    feature_dict[new_name] = {}
                feature_dict[new_name][gene] = value


#%%
common_genes = set(uniform_gene[0])
for genes in uniform_gene[1:]:
    common_genes &= set(genes)

common_feature = {}
for feature_name in feature_dict:
    common_feature[feature_name] = {}
    for gene in common_genes:
        common_feature[feature_name][gene] = feature_dict[feature_name][gene]

feature_matrix = np.zeros((len(common_genes), len(feature_dict.keys())))
for i,gene in enumerate(common_genes):
    for j,feature_name in enumerate(common_feature.keys()):
        feature_matrix[i, j] = common_feature[feature_name][gene]
tissue_expression_matrix = cal_tissue_matrix(expression_martrix, tissues_list,common_genes)
Weight = cal_nmf_weight(feature_matrix,tissue_expression_matrix)

weight_dict = {}
tissues = list(expression_martrix.keys())
feature = list(common_feature.keys())
for i, value in enumerate(Weight):
    weight_dict[tissues[i]] = {}
    for j, v in enumerate(Weight[i]):
        weight_dict[tissues[i]][feature[j]] = v
weight_dict = pd.DataFrame(weight_dict).T
weight_dict.to_csv(csv_file)
