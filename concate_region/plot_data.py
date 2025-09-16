#%%
import csv
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve, auc
import pandas as pd
import seaborn as sns
import argparse
import os

parser = argparse.ArgumentParser(description="Draw boxplot and ROC curve for two groups of data")
parser.add_argument("--group1_file", type=str, required=True, help="Path to group1 tsv file")
parser.add_argument("--group2_file", type=str, required=True, help="Path to group2 tsv file")
parser.add_argument("--group1_name", type=str, required=True, help="Name of group1 (e.g. cancer)")
parser.add_argument("--group2_name", type=str, required=True, help="Name of group2 (e.g. healthy)")
parser.add_argument("--title", type=str, default="Data", help="Title for the plots")
parser.add_argument("--save_dir", type=str, required=True, help="Save result directory")
args = parser.parse_args()

group1_file = args.group1_file
group2_file = args.group2_file
group1_name = args.group1_name
group2_name = args.group2_name
title = args.title
save_dir = args.save_dir


# group1_file = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/concate_region/pfe/cancer.pfe.tsv'
# group2_file = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/concate_region/pfe/healthy.pfe.tsv'
# group1_name = 'cancer'
# group2_name = 'healthy'
# title = ' pfe'


os.makedirs(save_dir, exist_ok=True)
group1_data = []
group2_data = []

f = open(group1_file,'r')
r = csv.reader(f,delimiter='\t')
for i in r:
    sample_name = i[0]
    value = float(i[1])
    group1_data.append(value)
f.close()

f = open(group2_file,'r')
r = csv.reader(f,delimiter='\t')
for i in r:
    sample_name = i[0]
    value = float(i[1])
    group2_data.append(value)
f.close()


# ==== 用 pandas+seaborn 整理数据 ====
df = pd.DataFrame({
    "Value": group1_data + group2_data,
    "Group": [group1_name] * len(group1_data) + [group2_name] * len(group2_data)
})

# ==== 画 seaborn 箱形图 ====
plt.figure(figsize=(6, 5))
sns.boxplot(x="Group", y="Value", data=df, palette="Set2")
sns.swarmplot(x="Group", y="Value", data=df, color=".25", size=3)  # 添加散点，更直观看分布
plt.title(f"Boxplot of {title}")
plt.tight_layout()
png_dir = os.path.join(save_dir, f'Boxplot.png')
plt.savefig(png_dir)


# ==== 画 ROC 曲线 ====
# 合并数据，构造标签：cancer=1, healthy=0
if np.mean(group1_data) > np.mean(group2_data):
    y_true = np.array([1] * len(group1_data) + [0] * len(group2_data))
    y_score = np.array(group1_data + group2_data)
else:
    y_true = np.array([0] * len(group1_data) + [1] * len(group2_data))
    y_score = np.array(group1_data + group2_data)

fpr, tpr, thresholds = roc_curve(y_true, y_score)
roc_auc = auc(fpr, tpr)

plt.figure(figsize=(6, 5))
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f"ROC curve (AUC = {roc_auc:.3f})")
plt.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title(f"ROC Curve of {title}")
plt.legend(loc="lower right")
plt.tight_layout()
png_dir = os.path.join(save_dir, f'ROC.png')
plt.savefig(png_dir)