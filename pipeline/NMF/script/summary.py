#%%
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


compare_weight_file = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/pos_healthy/nmf/nmf/compare.lst'
patient_weight_file = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/pos_healthy/nmf/nmf/patient.lst'
save_dir = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/pos_healthy/nmf/summary'


f = open(patient_weight_file)
patient_file_list = []
for i in f.readlines():
    patient_file_list.append(i.strip('\n'))


f = open(compare_weight_file)
compare_file_list = []
for i in f.readlines():
    compare_file_list.append(i.strip('\n'))


#%%
import pandas as pd
import os

tissue_list = set()
feature_list = set()


compare_data = {}
compare_list = set()
for its,file in enumerate(compare_file_list):
    file_name = os.path.basename(file)
    compare_list.add(file_name)
    print(file)
    data = pd.read_csv(file,index_col=0)
    compare_data[file_name] = {}
    for feature in data.columns:
        feature_list.add(feature)
        compare_data[file_name][feature] = {}
        for tissue in data.index:
            tissue_list.add(tissue)
            compare_data[file_name][feature][tissue] = data.loc[tissue][feature]


patient_data = {}
patient_list = set()
for its,file in enumerate(patient_file_list):
    file_name = os.path.basename(file)
    patient_list.add(file_name)
    print(file)
    data = pd.read_csv(file,index_col=0)
    patient_data[file_name] = {}
    for feature in data.columns:
        feature_list.add(feature)
        patient_data[file_name][feature] = {}
        for tissue in data.index:
            tissue_list.add(tissue)
            patient_data[file_name][feature][tissue] = data.loc[tissue][feature]
#%%
patient_summary_data = {}
for patient in patient_data:
    for tissue in tissue_list:
        if tissue not in patient_summary_data:
            patient_summary_data[tissue] = {}
        for feature in feature_list:
            if feature not in patient_summary_data[tissue]:
                patient_summary_data[tissue][feature] = 0

            temp = []
            for healthy in compare_data:
                temp.append(compare_data[healthy][feature][tissue])
            mean_value = np.mean(temp)
            std_value = np.std(temp)
            value = (patient_data[patient][feature][tissue]-mean_value)/std_value
            if value > 1:
                patient_summary_data[tissue][feature]+=1
            if patient_summary_data[tissue][feature] >= 7:
                print(patient,tissue, feature)
#%%
import pandas as pd
df = pd.DataFrame(patient_summary_data).T
import matplotlib.pyplot as plt
import seaborn as sns

df_filtered = df[df.max(axis=1) >= 5]

plt.figure(figsize=(12, 2))
sns.heatmap(
    df_filtered,
    cmap="Reds",          # 红色梯度，表示异常病人数量
    linewidths=0.5,
    linecolor="gray",
)

plt.xlabel("Feature")
plt.ylabel("Tissue")
plt.title("Patient Feature Enrichment Heatmap")
plt.tight_layout()
plt.show()

#%%
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
import numpy as np

X = []
y = []
feature = '19'
# patient = 1
for patient in patient_data:
    for tissue in df_filtered.index:
        X.append([
            patient_data[patient][feature][tissue]
            for feature in feature_list
        ])
        y.append(1)

# healthy = 0
feature = '19'
for compare in compare_data:
    for tissue in df_filtered.index:
        X.append([
            compare_data[compare][feature][tissue]
            for feature in feature_list
        ])
        y.append(0)
#%%
X = np.array(X)
y = np.array(y)


scaler = StandardScaler()
X = scaler.fit_transform(X)
clf = LogisticRegression(
    penalty='l2',
    solver='liblinear',
    class_weight='balanced',
    max_iter=1000
)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

y_score = cross_val_predict(
    clf, X, y,
    cv=cv,
    method='predict_proba'
)[:, 1]

auc = roc_auc_score(y, y_score)
print("ROC AUC:", auc)

from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt

# 计算 ROC 曲线
fpr, tpr, thresholds = roc_curve(y, y_score)

# 绘图
plt.figure(figsize=(6, 6))
plt.plot(fpr, tpr, label=f"ROC AUC = {auc:.3f}")
plt.plot([0, 1], [0, 1], linestyle='--', color='gray')

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve (Logistic Regression)")
plt.legend(loc="lower right")
plt.tight_layout()
plt.show()
