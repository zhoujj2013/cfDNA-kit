import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import argparse
import os
# 1. Create argument parser
parser = argparse.ArgumentParser(description="NMF multi-feature analysis script")
# 2. Add arguments
parser.add_argument("--sample_file", type=str, required=True,
                    help="Path to the sample list file, e.g., /path/to/samples.lst")
parser.add_argument("--save_dir", type=str, required=True,
                    help="Directory to save results")
parser.add_argument("--select_feature", type=str, required=True,
                    help="Feature to select (draw roc), e.g., fft_7")
parser.add_argument("--select_tissue", type=str, required=True,
                    help="Tissue to select (draw roc), e.g., Pancreas")

# 3. Parse arguments
args = parser.parse_args()

# 4. Use the arguments
sample_file = args.sample_file
save_dir = args.save_dir
select_feature = args.select_feature
select_tissue = args.select_tissue

# sample_file = '/mnt/dfc_data2/project/linyusen/project/31_cfdna_wps/project/NMF/multi_feature/samples.lst'
# save_dir = '/mnt/dfc_data2/project/linyusen/project/31_cfdna_wps/project/NMF/multi_feature/result'
# select_feature = 'fft_7'
# select_tissue = 'Pancreas'


os.makedirs(save_dir, exist_ok=True)
f = open(sample_file,'r')
sample_dict = {}
for i in f.readlines():
    i = i.strip('\n')
    group,name,path = i.split(' ')
    if group not in sample_dict:
        sample_dict[group] = {}
    sample_dict[group][name] = path

#%%
if 'control' not in sample_dict:
    raise KeyError("'control' key is missing in sample_dict")

group = 'control'
control_data = {}
tissue_list = set()
feature_list = set()
control_list = set()
for name in sample_dict[group]:
    control_list.add(name)
    file = sample_dict[group][name]
    data = pd.read_csv(file, index_col=0)
    control_data[name] = {}
    for feature in data.columns:
        feature_list.add(feature)
        control_data[name][feature] = {}
        for tissue in data.index:
            tissue_list.add(tissue)
            control_data[name][feature][tissue] = data.loc[tissue][feature]

tissue_list = list(tissue_list)
feature_list = list(feature_list)
control_list = list(control_list)
tissue_list.sort()
feature_list.sort()
control_list.sort()
experiment_data = {}
for group in sample_dict:
    if group == 'control':
        continue
    experiment_data[group] = {}
    for name in sample_dict[group]:
        file = sample_dict[group][name]
        data = pd.read_csv(file, index_col=0)
        experiment_data[group][name] = data

    patient_data = {}
    for name in experiment_data[group]:
        patient_weight = experiment_data[group][name]
        for tissue in tissue_list:
            if tissue not in patient_data:
                patient_data[tissue] = {}
            for feature in feature_list:
                if feature not in patient_data[tissue]:
                    patient_data[tissue][feature] = 0
                temp = []
                for patient in control_list:
                    temp.append(control_data[patient][feature][tissue])
                mean_value = np.mean(temp)
                std_value = np.std(temp)
                value = (patient_weight.loc[tissue][feature] - mean_value) / std_value
                if value > 1.5 or value < -1.50:
                    patient_data[tissue][feature] += 1
    control_df = pd.DataFrame(patient_data).T
    bwith = 3
    plt.figure(figsize=(12, 8))
    sns.heatmap(
        control_df,
        cmap="coolwarm",  # 蓝白红配色
        center=0,  # 中心对齐 0
        linewidths=0.5,
        linecolor="gray",
        annot=True,  # 显示数值
        fmt=".2f"  # 保留两位小数
    )
    plt.title("Heatmap")
    plt.ylabel("Tissue")
    plt.xlabel("Feature")
    plt.tick_params(axis='both', which='major', width=bwith, length=bwith * 3, labelsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'Feature_Zscore_Heatmap.png'))
    control_df.to_csv(os.path.join(save_dir, f'Feature_Zscore_Heatmap.csv'))



    experiment_value_dict = {}
    control_value_dict = {}
    for tissue in tissue_list:
        experiment_value_dict[tissue] = {}
        for name in experiment_data[group]:
            experiment_value_dict[tissue][name] = experiment_data[group][name].loc[tissue][select_feature]
    for tissue in tissue_list:
        control_value_dict[tissue] = {}
        for name in control_data:
            control_value_dict[tissue][name] = control_data[name][select_feature][tissue]
    experiment_df = pd.DataFrame(experiment_value_dict).T
    control_df = pd.DataFrame(control_value_dict).T
    experiment_long = experiment_df.reset_index().melt(id_vars='index', var_name='Patient', value_name='Value')
    experiment_long['Group'] = group
    control_long = control_df.reset_index().melt(id_vars='index', var_name='Patient', value_name='Value')
    control_long['Group'] = 'Control'
    combined_long = pd.concat([experiment_long, control_long], ignore_index=True)
    combined_long.rename(columns={'index': 'Tissue'}, inplace=True)
    plt.figure(figsize=(16, 8))
    sns.boxplot(data=combined_long, x='Tissue', y='Value', hue='Group', palette='Set2')
    bwith = 3
    ax = plt.gca()
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['bottom'].set_linewidth(bwith)
    ax.set_xlabel('Tissue', fontsize=16)
    ax.set_ylabel('Value', fontsize=16)
    plt.tick_params(axis='both', which='major', width=bwith, length=bwith * 3, labelsize=16)
    plt.xticks(rotation=90)
    plt.title(f'Expression Distribution per Tissue ({group} vs Control)', fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'boxplot.png'))
    combined_long.to_csv(os.path.join(save_dir, f'boxplot.csv'))

    combined_data = {}
    for tissue in control_value_dict:
        if tissue not in combined_data:
            combined_data[tissue] = {}
        for patient in control_value_dict[tissue]:
            combined_data[tissue][f'{group}_{patient}'] = control_value_dict[tissue][patient]
    for tissue in experiment_value_dict:
        for patient in experiment_value_dict[tissue]:
            combined_data[tissue][f'Control_{patient}'] = experiment_value_dict[tissue][patient]
    df = pd.DataFrame(combined_data).T
    df = df.dropna(axis=1, how="all")
    df = df.apply(lambda x: x.fillna(x.mean()), axis=0)
    mask_rows = df.std(axis=1) != 0
    df = df.loc[mask_rows, :]
    # 4. 仍然建议去掉标准差为 0 的列，避免 z_score 出 NaN
    mask_cols = df.std(axis=0) != 0
    df = df.loc[:, mask_cols.to_numpy()]
    col_colors = df.columns.to_series().map(lambda x: "red" if x.startswith(f"{group}_") else "blue")
    sns.clustermap(
        df,
        cmap="RdBu_r",
        figsize=(9, 8),
        z_score=0,
        col_colors=col_colors,
    )
    plt.savefig(os.path.join(save_dir, f'clustermap.png'))
    df.to_csv(os.path.join(save_dir, f'clustermap.csv'))



    scores1 = []
    for name in experiment_value_dict[select_tissue]:
        scores1.append(experiment_value_dict[select_tissue][name])
    scores2 = []
    for name in control_value_dict[select_tissue]:
        scores2.append(control_value_dict[select_tissue][name])
    scores = []
    labels = []
    if np.mean(scores1) > np.mean(scores2):
        for i in scores1:
            scores.append(i)
            labels.append(1)
        for i in scores2:
            scores.append(i)
            labels.append(0)
    else:
        for i in scores1:
            scores.append(i)
            labels.append(0)
        for i in scores2:
            scores.append(i)
            labels.append(1)

    scores = np.array(scores)
    labels = np.array(labels)
    # 2. 遍历阈值找到最优 Youden's J
    thresholds = np.unique(scores)
    best_threshold = None
    best_score = -np.inf

    for th in thresholds:
        preds = (scores > th).astype(int)
        TP = np.sum((preds == 1) & (labels == 1))
        FP = np.sum((preds == 1) & (labels == 0))
        TN = np.sum((preds == 0) & (labels == 0))
        FN = np.sum((preds == 0) & (labels == 1))

        sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
        specificity = TN / (TN + FP) if (TN + FP) > 0 else 0

        J = sensitivity + specificity - 1

        if J > best_score:
            best_score = J
            best_threshold = th

    # 3. 计算最终混淆矩阵
    preds = (scores > best_threshold).astype(int)
    TP = np.sum((preds == 1) & (labels == 1))
    FP = np.sum((preds == 1) & (labels == 0))
    TN = np.sum((preds == 0) & (labels == 0))
    FN = np.sum((preds == 0) & (labels == 1))

    accuracy = (TP + TN) / (TP + TN + FP + FN)
    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    f1_score = (2 * TP) / (2 * TP + FP + FN) if (TP + FP + FN) > 0 else 0

    # ===== 图 1：混淆矩阵 =====
    plt.figure(figsize=(8, 8))
    matrix = np.array([[TP, FN],
                       [FP, TN]])

    sns.heatmap(matrix, annot=True, fmt='d', cmap='Blues', cbar=False,
                xticklabels=['Predicted Patient', 'Predicted Healthy'],
                yticklabels=['Actual Patient', 'Actual Healthy'],
                annot_kws={"size": 18})

    plt.title(f'Confusion Matrix (Threshold={best_threshold:.3f})', fontsize=14)
    plt.ylabel('Actual')
    plt.xlabel('Predicted')
    matrix_df = pd.DataFrame(matrix,
                             index=['Actual Patient', 'Actual Healthy'],
                             columns=['Predicted Patient', 'Predicted Healthy'])
    matrix_df .to_csv(os.path.join(save_dir, f'matrix.csv'))

    # 在热图下方显示指标
    metrics_text = (
        f"Accuracy: {accuracy:.3f}    "
        f"Sensitivity: {sensitivity:.3f}    "
        f"Specificity: {specificity:.3f}\n"
        f"Precision: {precision:.3f}    "
        f"F1 Score: {f1_score:.3f}    "
        f"Youden's J: {best_score:.3f}"
    )
    plt.text(0.5, -0.15, metrics_text, ha='center', va='top', fontsize=12, transform=plt.gca().transAxes)

    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'matrix.png'))

    # ===== 图 2：ROC 曲线 =====
    plt.figure(figsize=(8, 8))
    fpr, tpr, roc_thresholds = roc_curve(labels, scores)
    roc_auc = auc(fpr, tpr)

    plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC curve (AUC = {roc_auc:.3f})')
    plt.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=16)
    plt.ylabel('True Positive Rate', fontsize=16)
    plt.title('ROC Curve', fontsize=14)
    plt.legend(loc="lower right")

    # 设置边框
    bwith = 3
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(bwith)

    # 设置刻度
    plt.tick_params(axis='both', which='major', width=bwith, length=bwith * 3, labelsize=16)

    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'ROC.png'))
