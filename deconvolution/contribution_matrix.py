#%%
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Process NMF weights and save normalized contributions")
parser.add_argument(
    "--compare_weight_lst",
    type=str,
    required=True,
    help="Path to the compare.lts file"
)
parser.add_argument(
    "--patient_weight",
    type=str,
    required=True,
    help="Path to the patient weights CSV file"
)
parser.add_argument(
    "--save_csv",
    type=str,
    required=True,
    help="Path to save the normalized contribution results"
)
args = parser.parse_args()
compare_weight_lst = args.compare_weight_lst
patient_weight = args.patient_weight
save_csv = args.save_csv

# compare_weight_lst = '/mnt/dfc_data1/home/linyusen/31_cfdna_wps/project/NMF/compare.lts'
# patient_weight = f'/mnt/dfc_data2/project/linyusen/database/46_cfdna/nmf_data/PL230613012weight.csv'
# save_csv = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/nmf_data/normalization.contribution.csv'

f = open(compare_weight_lst)
compare_file_list = []
for i in f.readlines():
    compare_file_list.append(i[:-1])

compare_data = {}

tissue_list = set()
feature_list = set()
compare_list = set()
for its,file in enumerate(compare_file_list):
    compare_list.add(its)
    print(file)
    data = pd.read_csv(file,index_col=0)
    compare_data[its] = {}
    for feature in data.columns:
        feature_list.add(feature)
        compare_data[its][feature] = {}
        for tissue in data.index:
            tissue_list.add(tissue)
            compare_data[its][feature][tissue] = data.loc[tissue][feature]
#%%
patient_weight = pd.read_csv(patient_weight,index_col=0)
patient_data = {}
for tissue in tissue_list:
    patient_data[tissue] = {}
    for feature in feature_list:
        temp = []
        for patient in compare_list:
            temp.append(compare_data[patient][feature][tissue])
        mean_value = np.mean(temp)
        std_value = np.std(temp)
        value = (patient_weight.loc[tissue][feature]-mean_value)/std_value
        patient_data[tissue][feature] = value
#%%
# 转换字典为 DataFrame
cancer_df = pd.DataFrame(patient_data).T
cancer_df.to_csv(save_csv)