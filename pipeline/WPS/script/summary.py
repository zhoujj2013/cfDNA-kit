import os
import pickle
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import argparse

parse = argparse.ArgumentParser()
parse.add_argument('-single_cell_matrix',required=True)
parse.add_argument('-data_path',required=True)
parse.add_argument('-save_path',required=True)
args = parse.parse_args()
single_cell_matrix = args.single_cell_matrix
data_path = args.data_path
save_path = args.save_path

# single_cell_matrix = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/nerous_all_cell/single_cell/single_cell.expression.csv'
# data_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/data2/'
# save_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/nerous_all_cell/data'

df = pd.read_csv(single_cell_matrix,header=0,index_col=0)
dict_data = df.to_dict('dict')
#%%
gene_intensity_198 = {}
gene_intensity_ratio = {}
if os.path.exists(save_path) == False:
    os.makedirs(save_path)
for patient in os.listdir(data_path):
    print(patient)
    gene_intensity_198[patient] = {}
    gene_intensity_ratio[patient] = {}
    fft_file = os.path.join(data_path,patient,'fft.pickle')
    with open(fft_file, 'rb') as file:
        fft = pickle.load(file)
        for gene in fft.keys():
            temp_198 = []
            for fre in range(190,200):
                if fre in fft[gene].keys():
                    temp_198.append(fft[gene][fre])
            temp_225 = []
            for fre in range(215,225):
                if fre in fft[gene].keys():
                    temp_225.append(fft[gene][fre])
            if len(temp_198) >= 1:
                gene_intensity_198[patient][gene] = np.mean(temp_198)
            if len(temp_198) >= 1 and len(temp_225) >= 1:
                gene_intensity_ratio[patient][gene] = np.mean(temp_198) / np.mean(temp_225)
# 将字典保存到 pickle 文件
with open(os.path.join(save_path,'gene_intensity_ratio.pkl'), 'wb') as pickle_file:
    pickle.dump(gene_intensity_ratio, pickle_file)
with open(os.path.join(save_path,'gene_intensity_198.pkl'), 'wb') as pickle_file:
    pickle.dump(gene_intensity_198, pickle_file)
#%%
patient_cell_corelation_198 = {}
for patient in gene_intensity_198.keys():
    print(patient)
    cell_corelation = {}
    for cell in dict_data.keys():
        x1 = []
        x2 = []
        gene_list = list(set(gene_intensity_198[patient].keys())&set(dict_data[cell].keys()))
        for gene in gene_list:
            x1.append(gene_intensity_198[patient][gene])
            x2.append(dict_data[cell][gene])
        correlation_coefficient, p_value = pearsonr(x1, x2)
        cell_corelation[cell] = correlation_coefficient
    patient_cell_corelation_198[patient] = cell_corelation
with open(os.path.join(save_path,'patient_cell_corelation_198.pkl'),'wb') as pickle_file:
    pickle.dump(patient_cell_corelation_198, pickle_file)
#%%
patient_cell_corelation_ratio = {}
for patient in gene_intensity_ratio.keys():
    print(patient)
    cell_corelation = {}
    for cell in dict_data.keys():
        x1 = []
        x2 = []
        gene_list = list(set(gene_intensity_ratio[patient].keys())&set(dict_data[cell].keys()))
        for gene in gene_list:
            x1.append(gene_intensity_ratio[patient][gene])
            x2.append(dict_data[cell][gene])
        correlation_coefficient, p_value = pearsonr(x1, x2)
        cell_corelation[cell] = correlation_coefficient
    patient_cell_corelation_ratio[patient] = cell_corelation
with open(os.path.join(save_path,'patient_cell_corelation_ratio.pkl'),'wb') as pickle_file:
    pickle.dump(patient_cell_corelation_ratio, pickle_file)

