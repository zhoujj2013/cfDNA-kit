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
    help="Path to the pancreas_control.lst file"
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
patient_weight_path= patient_weight
# compare_weight_lst = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/pos_healthy/nmf/nmf/compare.lst'
# patient_weight = f'/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/pos_healthy/nmf/nmf/pso8.weight.csv'
# save_csv = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/pos_healthy/nmf/contribution/pso8.normalization.contribution.csv'

f = open(compare_weight_lst)
compare_file_list = []
for i in f.readlines():
    compare_file_list.append(i.strip('\n'))
if patient_weight in compare_file_list:
    compare_file_list.remove(patient_weight)

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

tissue_list = list(tissue_list)
feature_list = list(feature_list)
compare_list = list(compare_list)

tissue_list.sort()
feature_list.sort()
compare_list.sort()

#%%
patient_weight = pd.read_csv(patient_weight,index_col=0)
patient_data = {}
for tissue in tissue_list:
    # _,cell,_ = tissue.split('#')
    # if cell in ['CD4-positive-helper-T-cell','memory-B-cell','CD141-positive-myeloid-dendritic-cell',
    #             'naive-B-cell','Langerhans-cell','naive-thymus-derived-CD4-positive-alpha-beta-T-cell']:
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

import seaborn as sns
import matplotlib.pyplot as plt

# 如果希望所有特征的颜色在同一 scale 下，可以计算全局的 vmin 和 vmax
vmin = -1.96
vmax = 1.96

plt.figure(figsize=(12, 8))
sns.heatmap(
    cancer_df,
    cmap="coolwarm",      # 蓝白红配色
    center=0,             # 中心对齐 0
    linewidths=0.5,
    linecolor="gray",
    vmin=vmin,
    vmax=vmax,
    annot=False,           # 显示数值
    # fmt=".2f"             # 保留两位小数
)
import os
patient_name = os.path.basename(patient_weight_path)
plt.title(patient_name)
plt.ylabel("Tissue")
plt.xlabel("Feature")
plt.tight_layout()
plt.savefig(save_csv+'.png')