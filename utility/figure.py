import os.path
import pickle
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression
from scipy import stats
import seaborn as sns
from sklearn.manifold import TSNE
import csv
import pandas as pd
import argparse

parse = argparse.ArgumentParser()
parse.add_argument('-gene_intensity_path',required=True)
parse.add_argument('-patient_cell_corelation_path',required=True)
parse.add_argument('-single_cell_expression',required=True)
parse.add_argument('-sample_file',required=True)
parse.add_argument('-save_path',required=True)
args = parse.parse_args()
gene_intensity_ratio_path = args.gene_intensity_ratio_path
patient_cell_corelation_ratio_path = args.patient_cell_corelation_ratio_path
single_cell_expression = args.single_cell_expression
sample_file = args.sample_file
save_path = args.save_path


# gene_intensity_ratio_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/pso_hc/summary/gene_intensity_ratio.pkl'
# patient_cell_corelation_ratio_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/pso_hc/summary/patient_cell_corelation_ratio.pkl'
# single_cell_expression = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/pso_hc/single_cell/single_cell.expression.csv'
# sample_file = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/pso_hc/sample.lst'
# save_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/pso_hc/figure'
#%%
if os.path.exists(save_path) == False:
    os.makedirs(save_path)
def min_max_normalization(data):
    min_val = np.min(data)
    max_val = np.max(data)
    normalized_data = (data - min_val) / (max_val - min_val)
    return normalized_data
def calculate_variance(data):
    mean = np.mean(data)
    variance = np.mean((data - mean) ** 2)
    return variance
with open(gene_intensity_ratio_path, 'rb') as file:
    gene_intensity = pickle.load(file)
with open(patient_cell_corelation_ratio_path, 'rb') as file:
    patient_cell_corelation = pickle.load(file)
for patient in patient_cell_corelation:
    for cell in patient_cell_corelation[patient]:
        patient_cell_corelation[patient][cell] = abs(patient_cell_corelation[patient][cell])
sample_list = {}
with open(sample_file,'r') as f:
    reader = csv.reader(f,delimiter='\t')
    for line in reader:
        type,group,name,path = line
        if type not in sample_list:
            sample_list[type] = {}
        sample_list[type][name] = {'path':path,'group':group}
control_group = []
control_group_name= ''
for type in sample_list.keys():
    for patient in sample_list[type].keys():
        path = sample_list[type][patient]['path']
        group = sample_list[type][patient]['group']
        if group == 'control':
            control_group.append(patient)
            control_group_name = type

col_gene = gene_intensity[list(gene_intensity.keys())[0]].keys()
for patient in gene_intensity.keys():
    col_gene = set(col_gene) & set(gene_intensity[patient].keys())
gene_intensity_sort = {}
for sample in gene_intensity.keys():
    # 对字典按值排序（升序）
    sorted_items = sorted(gene_intensity[sample].items(), key=lambda x: x[1],reverse=True)
    # 创建一个新字典，记录每个键在排序后的序号
    rank_dict = {k: idx for idx, (k, v) in enumerate(sorted_items, 1)}
    gene_intensity_sort[sample] = rank_dict

#%%
feature = {}
activate_gene = {}
deactivate_gene = {}
comapre_group_name = list(sample_list.keys())
comapre_group_name.remove(control_group_name)
for type in comapre_group_name:
    print(type)
    activate_gene[type] = []
    deactivate_gene[type] = []
    compare_group = list(sample_list[type].keys())
    for index,gene in enumerate(col_gene):
        compare = []
        control = []

        compare_sort = []
        control_sort = []
        for patient in compare_group:
            compare.append(gene_intensity[patient][gene])
        for patient in control_group:
            control.append(gene_intensity[patient][gene])

        for patient in compare_group:
            compare_sort.append(gene_intensity_sort[patient][gene])
        for patient in control_group:
            control_sort.append(gene_intensity_sort[patient][gene])

        fc = np.log2(np.mean(control) / np.mean(compare))
        var = np.var(compare)
        temp = {
            'gene':gene,
            'var': var,
            'fc': fc,
            'fft_hc':np.mean(control),
            'fft_ds': np.mean(compare),
            'fft_hc_s': np.mean(control_sort),
            'fft_ds_s': np.mean(compare_sort)
        }
        # if var > 50:
        #     continue
        if abs(fc) > 0:
            if fc > 0:

                activate_gene[type].append(temp)
            elif fc < -0:
                deactivate_gene[type].append(temp)
            for ii in gene_intensity.keys():
                if ii not in feature:
                    feature[ii] = {}
                if gene not in feature[ii].keys():
                    feature[ii][gene] = gene_intensity_sort[ii][gene]
de_activate_gene = {}
zscore_gene = {}
comapre_group_name = list(sample_list.keys())
comapre_group_name.remove(control_group_name)
for type in comapre_group_name:
    print(type)
    de_activate_gene[type] = []
    compare_group = list(sample_list[type].keys())
    for index,gene in enumerate(col_gene):
        compare = []
        control = []

        compare_sort = []
        control_sort = []
        for patient in compare_group:
            compare.append(gene_intensity[patient][gene])
        for patient in control_group:
            control.append(gene_intensity[patient][gene])

        for patient in compare_group:
            if patient not in zscore_gene:
                zscore_gene[patient] = {}
            compare_sort.append(gene_intensity_sort[patient][gene])
        for patient in control_group:
            if patient not in zscore_gene:
                zscore_gene[patient] = {}
            control_sort.append(gene_intensity_sort[patient][gene])

        meam_control_sort = np.mean(control_sort)
        std_control_sort = np.std(control_sort)
        for patient in compare_group:
            zscore_gene[patient][gene] = (gene_intensity_sort[patient][gene] - meam_control_sort)/std_control_sort
        for patient in control_group:
            zscore_gene[patient][gene] = (gene_intensity_sort[patient][gene] - meam_control_sort)/std_control_sort
        fc = np.log2(np.mean(control) / np.mean(compare))
        var = np.var(compare)
        temp = {
            'gene':gene,
            'var': var,
            'fc': fc,
            'fft_hc':np.mean(control),
            'fft_ds': np.mean(compare),
            'fft_hc_s': np.mean(control_sort),
            'fft_ds_s': np.mean(compare_sort)
        }
        # if var > 50:
        #     continue

        de_activate_gene[type].append(temp)
#%%
f = open('/mnt/dfc_data2/project/linyusen/database/02_hg/hg38/ensg_offic_gene')
human_ensg_gene = {}
for i in f.readlines():
    ensg,gene = i[:-1].split('\t')
    human_ensg_gene[ensg] = gene

for type in de_activate_gene.keys():
    activate_save_path = os.path.join(save_path,type+f'.gene.txt')
    print(activate_save_path)

    patient_list = []
    for patient in gene_intensity:
        if type in patient:
            patient_list.append(patient)
    patient_sort_list = []
    for patient in gene_intensity:
        if type in patient:
            patient_sort_list.append(patient+'_sort')
    patient_zscore_list = []
    for patient in gene_intensity:
        if type in patient:
            patient_zscore_list.append(patient + '_zscore')

    healthy_list = []
    for patient in control_group:
        healthy_list.append(patient)
    healthy_sort_list = []
    for patient in control_group:
        healthy_sort_list.append(patient + '_sort')
    healthy_zscore_list = []
    for patient in control_group:
        healthy_zscore_list.append(patient + '_zscore')
    f = open(activate_save_path, 'w')
    aaa = '\t'.join(patient_list)
    bbb = '\t'.join(patient_sort_list)
    ccc = '\t'.join(healthy_list)
    ddd = '\t'.join(healthy_sort_list)
    eee = '\t'.join(patient_zscore_list)
    fff = '\t'.join(healthy_zscore_list)
    f.write(f'ensg\tgene\tfc\tvar\tfft_healthy\tfft_disease\tfft_sort_healthy\tfft_sort_disease\t{aaa}\t{bbb}\t{ccc}\t{ddd}\t{eee}\t{fff}\n')
    for temp in de_activate_gene[type]:
        ensg = temp['gene']
        if ensg in human_ensg_gene:
            gene = human_ensg_gene[ensg]
        else:
            gene = ensg
        fc = temp['fc']
        var = temp['var']
        fft_hc = temp['fft_hc']
        fft_ds = temp['fft_ds']
        fft_hc_s = temp['fft_hc_s']
        fft_ds_s = temp['fft_ds_s']
        fft_patient = []
        for patient in patient_list:
            fft_patient.append(str(gene_intensity[patient][ensg]))
        fft_patient_sort = []
        for patient in patient_list:
            fft_patient_sort.append(str(gene_intensity_sort[patient][ensg]))
        aaa = '\t'.join(fft_patient)
        bbb = '\t'.join(fft_patient_sort)

        fft_healthy = []
        for patient in control_group:
            fft_healthy.append(str(gene_intensity[patient][ensg]))
        fft_healthy_sort = []
        for patient in control_group:
            fft_healthy_sort.append(str(gene_intensity_sort[patient][ensg]))
        ccc = '\t'.join(fft_healthy)
        ddd = '\t'.join(fft_healthy_sort)

        fft_patient_zscore = []
        for patient in patient_list:
            fft_patient_zscore.append(str(zscore_gene[patient][ensg]))
        fft_healthy_zscore = []
        for patient in control_group:
            fft_healthy_zscore.append(str(zscore_gene[patient][ensg]))
        eee = '\t'.join(fft_healthy)
        fff = '\t'.join(fft_healthy_zscore)

        f.write(f'{ensg}\t{gene}\t{fc}\t{var}\t{fft_hc}\t{fft_ds}\t{fft_hc_s}\t{fft_ds_s}\t{aaa}\t{bbb}\t{ccc}\t{ddd}\t{eee}\t{fff}\n')
    f.close()
#%%
x = []
y = []
name_list = []
for index,type in enumerate(sample_list):
    for key in sample_list[type]:
        x.append(np.array(list(feature[key].values())))
        y.append(index)
        name_list.append(key)
for index,type in enumerate(sample_list):
    for key in sample_list[type]:
        x.append(np.array(list(feature[key].values())))
        y.append(index)
x = np.array(x)
num_colors = len(sample_list.keys())
cmap = plt.get_cmap('tab10')
colors = [cmap(i) for i in np.linspace(0, 1, num_colors)]
pls = PLSRegression(n_components=3, max_iter=1000)
X_train_pls = pls.fit_transform(x, y)[0]
fig = plt.figure(figsize=(10, 10), dpi=500)
ax = fig.add_subplot(111, projection='3d')
for label, color in zip(np.unique(y), colors):
    ax.scatter(X_train_pls[y == label, 0],
               X_train_pls[y == label, 1],
               X_train_pls[y == label, 2],
               s=100,
               label=f'{list(sample_list.keys())[label]}',
               color=color,
               alpha=0.8,
               marker='x',
               linewidths=4)
bwith = 3
plt.tick_params(axis='both', which='major', width=bwith, length=bwith*3, labelsize=16)
ax.set_title('PLS-DA of Dataset')
ax.legend(bbox_to_anchor=(0, 1), loc='upper right')
plt.savefig(os.path.join(save_path,'cluster_PLS-DA.png'))

x = []
y = []
name_list = []
for index,type in enumerate(sample_list):
    for key in sample_list[type]:
        x.append(np.array(list(feature[key].values())))
        y.append(index)
        name_list.append(key)

for index,type in enumerate(sample_list):
    for key in sample_list[type]:
        x.append(np.array(list(feature[key].values())))
        y.append(index)
x = np.array(x)
num_colors = len(sample_list.keys())
cmap = plt.get_cmap('tab10')
colors = [cmap(i) for i in np.linspace(0, 1, num_colors)]
pca = PCA(n_components=3)
X_pca = pca.fit_transform(x)
fig = plt.figure(figsize=(10, 10), dpi=500)
ax = fig.add_subplot(111, projection='3d')
for label, color in zip(np.unique(y), colors):
    ax.scatter(X_pca[y == label, 0],
               X_pca[y == label, 1],
               X_pca[y == label, 2],
               s=100,
               label=f'{list(sample_list.keys())[label]}',
               color=color,
               alpha=0.8,
               marker='x',
               linewidths=4)
bwith = 3
plt.tick_params(axis='both', which='major', width=bwith, length=bwith*3, labelsize=16)
ax.set_title('PCA of Dataset')
ax.legend(bbox_to_anchor=(0, 1), loc='upper right')
plt.savefig(os.path.join(save_path,'cluster_PCA.png'))

x = []
y = []
name_list = []
for index,type in enumerate(sample_list):
    for key in sample_list[type]:
        x.append(np.array(list(feature[key].values())))
        y.append(index)
        name_list.append(key)

for index,type in enumerate(sample_list):
    for key in sample_list[type]:
        x.append(np.array(list(feature[key].values())))
        y.append(index)
x = np.array(x)
num_colors = len(sample_list.keys())
cmap = plt.get_cmap('tab10')
colors = [cmap(i) for i in np.linspace(0, 1, num_colors)]
n_samples = x.shape[0]
perplexity = min(30, n_samples - 1)
tsne = TSNE(n_components=3, perplexity=perplexity, random_state=42)
X_tsne = tsne.fit_transform(x)
fig = plt.figure(figsize=(10, 10), dpi=500)
ax = fig.add_subplot(111, projection='3d')
for label, color in zip(np.unique(y), colors):
    ax.scatter(X_tsne[y == label, 0],
               X_tsne[y == label, 1],
               X_tsne[y == label, 2],
               s=100,
               label=f'{list(sample_list.keys())[label]}',
               color=color,
               alpha=0.8,
               marker='x',
               linewidths=4)
bwith = 3
plt.tick_params(axis='both', which='major', width=bwith, length=bwith*3, labelsize=16)
ax.set_title('TSNE of Dataset')
ax.legend(bbox_to_anchor=(0, 1), loc='upper right')
plt.savefig(os.path.join(save_path,'cluster_TSEN.png'))

#%%

x = []
y = []
name_list = []
for index,type in enumerate(sample_list):
    for key in sample_list[type]:
        x.append(np.array(list(feature[key].values())))
        y.append(index)
        name_list.append(key)

for index,type in enumerate(sample_list):
    for key in sample_list[type]:
        x.append(np.array(list(feature[key].values())))
        y.append(index)
x = np.array(x)
num_colors = len(sample_list.keys())
cmap = plt.get_cmap('tab10')
colors = [cmap(i) for i in np.linspace(0, 1, num_colors)]
pls = PLSRegression(n_components=3, max_iter=1000)
X_train_pls = pls.fit_transform(x, y)[0]
fig = plt.figure(figsize=(10, 10), dpi=500)
ax = fig.add_subplot(111, projection='3d')
for label, color in zip(np.unique(y), colors):
    ax.scatter(X_train_pls[y == label, 0],
               X_train_pls[y == label, 1],
               X_train_pls[y == label, 2],
               s=100,
               label=f'{list(sample_list.keys())[label]}',
               color=color,
               alpha=0.8,
               marker='x',
               linewidths=4)
for (x,y,z),(name) in zip(X_train_pls,name_list):
    ax.text(x,y,z,name)
bwith = 3
plt.tick_params(axis='both', which='major', width=bwith, length=bwith*3, labelsize=16)
ax.set_title('PLS-DA of Dataset')
ax.legend(bbox_to_anchor=(0, 1), loc='upper right')
plt.savefig(os.path.join(save_path,'cluster_PLS-DA.label.png'))


x = []
y = []
name_list = []
for index,type in enumerate(sample_list):
    for key in sample_list[type]:
        x.append(np.array(list(feature[key].values())))
        y.append(index)
        name_list.append(key)

for index,type in enumerate(sample_list):
    for key in sample_list[type]:
        x.append(np.array(list(feature[key].values())))
        y.append(index)
x = np.array(x)
num_colors = len(sample_list.keys())
cmap = plt.get_cmap('tab10')
colors = [cmap(i) for i in np.linspace(0, 1, num_colors)]
pca = PCA(n_components=3)
X_pca = pca.fit_transform(x)
fig = plt.figure(figsize=(10, 10), dpi=500)
ax = fig.add_subplot(111, projection='3d')
for label, color in zip(np.unique(y), colors):
    ax.scatter(X_pca[y == label, 0],
               X_pca[y == label, 1],
               X_pca[y == label, 2],
               s=100,
               label=f'{list(sample_list.keys())[label]}',
               color=color,
               alpha=0.8,
               marker='x',
               linewidths=4)
for (x,y,z),(name) in zip(X_train_pls,name_list):
    ax.text(x,y,z,name)
bwith = 3
plt.tick_params(axis='both', which='major', width=bwith, length=bwith*3, labelsize=16)
ax.set_title('PCA of Dataset')
ax.legend(bbox_to_anchor=(0, 1), loc='upper right')
plt.savefig(os.path.join(save_path,'cluster_PCA.label.png'))


x = []
y = []
name_list = []
for index,type in enumerate(sample_list):
    for key in sample_list[type]:
        x.append(np.array(list(feature[key].values())))
        y.append(index)
        name_list.append(key)

for index,type in enumerate(sample_list):
    for key in sample_list[type]:
        x.append(np.array(list(feature[key].values())))
        y.append(index)
x = np.array(x)
num_colors = len(sample_list.keys())
cmap = plt.get_cmap('tab10')
colors = [cmap(i) for i in np.linspace(0, 1, num_colors)]
n_samples = x.shape[0]
perplexity = min(30, n_samples - 1)
tsne = TSNE(n_components=3, perplexity=perplexity, random_state=42)
X_tsne = tsne.fit_transform(x)
fig = plt.figure(figsize=(10, 10), dpi=500)
ax = fig.add_subplot(111, projection='3d')
for label, color in zip(np.unique(y), colors):
    ax.scatter(X_tsne[y == label, 0],
               X_tsne[y == label, 1],
               X_tsne[y == label, 2],
               s=100,
               label=f'{list(sample_list.keys())[label]}',
               color=color,
               alpha=0.8,
               marker='x',
               linewidths=4)
for (x,y,z),(name) in zip(X_train_pls,name_list):
    ax.text(x,y,z,name)
bwith = 3
plt.tick_params(axis='both', which='major', width=bwith, length=bwith*3, labelsize=16)
ax.set_title('TSNE of Dataset')
ax.legend(bbox_to_anchor=(0, 1), loc='upper right')
plt.savefig(os.path.join(save_path,'cluster_TSEN.label.png'))


#%%
# 读取TSV文件
df = pd.read_csv(single_cell_expression,header=0,index_col=0)
dict_data = df.to_dict('dict')
tissue_correction = {}
for type in sample_list:
    tissue_correction[type] = {}
    for patient in sample_list[type]:
        tissue_correction_temp = {}
        for i in patient_cell_corelation[patient]:
            if i not in tissue_correction_temp:
                tissue_correction_temp[i] = []
            tissue_correction_temp[i].append(patient_cell_corelation[patient][i])
        tissue_correction_temp = dict(sorted(tissue_correction_temp.items(), key=lambda x: x[1],reverse=True))
        for index,i in enumerate(tissue_correction_temp):

            if i in dict_data.keys():
                celltype,cell,tissue = i.split('#')
                if tissue not in tissue_correction[type]:
                    tissue_correction[type][tissue] = []
                tissue_correction[type][tissue].append(index)
tissue_rank_dict = {}
for type in tissue_correction.keys():
    tissue_rank_dict[type] = dict(sorted(tissue_correction[type].items(), key=lambda x: sum(x[1]) / len(x[1])))

x_ticks = list(tissue_rank_dict[list(sample_list.keys())[0]].keys())
data = {}
for patient in sample_list.keys():
    data[patient] = {}
    for tissue in x_ticks:
        if patient not in data[patient]:
            data[patient][tissue] = []
        data[patient][tissue] = tissue_rank_dict[patient][tissue]
plot_data = []
group_labels = []
subgroup_labels = []
for group, subgroups in data.items():
    for subgroup, values in subgroups.items():
        plot_data.extend(values)
        group_labels.extend([group] * len(values))
        subgroup_labels.extend([subgroup] * len(values))
plt.figure(figsize=(12, 7), dpi=500)
sns.boxplot(x=subgroup_labels, y=plot_data, hue=group_labels)
# 设置图形属性
bwith = 3
ax = plt.gca()  # 获取边框
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
ax.spines['bottom'].set_linewidth(bwith)
# 移动y轴刻度到右边
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
plt.tick_params(axis='both', which='major', width=bwith, length=bwith*3, labelsize=16)
plt.tick_params(axis='x', rotation=90)
plt.tick_params(axis='y', rotation=90)
plt.ylabel('Ranking', size=16)
plt.title('tissue rank')
plt.legend(title='Dataset', bbox_to_anchor=(0, 1), loc='upper right')
plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.35)
plt.savefig(os.path.join(save_path,'tissue_rank.png'))
tissue_correction = {}
for type in sample_list.keys():
    if type not in tissue_correction:
        tissue_correction[type] = {}
    for patient in sample_list[type].keys():
        tissue_correction_temp = {}
        for cell_tissue in patient_cell_corelation[patient]:
            if cell_tissue not in tissue_correction_temp:
                tissue_correction_temp[cell_tissue] = []
            tissue_correction_temp[cell_tissue].append(patient_cell_corelation[patient][cell_tissue])
        tissue_correction_temp = dict(sorted(tissue_correction_temp.items(), key=lambda x: x[1],reverse=True))
        for index,cell_tissue in enumerate(tissue_correction_temp):
            celltype,cell, tissue = cell_tissue.split('#')
            if celltype != 'none' and celltype != 'None' and celltype != 'NONE':
                if celltype not in tissue_correction[type]:
                    tissue_correction[type][celltype] = []
                tissue_correction[type][celltype].append(index)
tissue_rank_dict = {}
for type in tissue_correction.keys():
    tissue_rank_dict[type] = dict(sorted(tissue_correction[type].items(), key=lambda x: sum(x[1]) / len(x[1])))
# 创建示例数据
x_ticks = list(tissue_rank_dict[list(sample_list.keys())[0]].keys())
data = {}
for patient in sample_list.keys():
    data[patient] = {}
    for tissue in x_ticks:
        if patient not in data[patient]:
            data[patient][tissue] = []
        data[patient][tissue] = tissue_rank_dict[patient][tissue]
plot_data = []
group_labels = []
subgroup_labels = []
for group, subgroups in data.items():
    for subgroup, values in subgroups.items():
        plot_data.extend(values)
        group_labels.extend([group] * len(values))
        subgroup_labels.extend([subgroup] * len(values))
plt.figure(figsize=(5, 7), dpi=500)
sns.boxplot(x=subgroup_labels, y=plot_data, hue=group_labels)
# 设置图形属性
bwith = 3
ax = plt.gca()  # 获取边框
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
ax.spines['bottom'].set_linewidth(bwith)
# 移动y轴刻度到右边
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
plt.tick_params(axis='both', which='major', width=bwith, length=bwith*3, labelsize=16)
plt.tick_params(axis='x', rotation=90)
plt.tick_params(axis='y', rotation=90)
plt.ylabel('Ranking', size=16)
plt.title('tissue rank')
plt.legend(title='Dataset', bbox_to_anchor=(0, 1), loc='upper right')
plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.35)
plt.savefig(os.path.join(save_path,'cell_type_rank.png'))
#%%
cell_tissue_group = {}
patient = list(patient_cell_corelation.keys())[0]
for cell_tissue in patient_cell_corelation[patient].keys():
    celltype,cell,tissue = cell_tissue.split('#')
    if cell not in cell_tissue_group.keys():
        cell_tissue_group[cell] = []
    cell_tissue_group[cell].append(cell_tissue)
for cell in list(cell_tissue_group.keys()):
    if len(cell_tissue_group[cell]) < 6:
        cell_tissue_group.pop(cell)
for cell in list(cell_tissue_group.keys()):
    tissue_correction = {}
    for patient_type in sample_list:
        tissue_correction[patient_type] = {}
        for patient in sample_list[patient_type]:
            tissue_correction_temp = {}
            for i in patient_cell_corelation[patient]:
                if i not in tissue_correction_temp:
                    tissue_correction_temp[i] = []
                tissue_correction_temp[i].append(patient_cell_corelation[patient][i])
            tissue_correction_temp = dict(sorted(tissue_correction_temp.items(), key=lambda x: x[1], reverse=True))
            for index, cell_tissue in enumerate(tissue_correction_temp):
                if cell_tissue in cell_tissue_group[cell]:
                    celltype, cell, tissue = cell_tissue.split('#')
                    if tissue not in tissue_correction[patient_type]:
                        tissue_correction[patient_type][tissue] = []
                    tissue_correction[patient_type][tissue].append(index)
    tissue_rank_dict = {}
    for type in tissue_correction.keys():
        tissue_rank_dict[type] = dict(sorted(tissue_correction[type].items(), key=lambda x: sum(x[1]) / len(x[1])))

    # 创建示例数据
    x_ticks = list(tissue_rank_dict[list(sample_list.keys())[0]].keys())
    data = {}
    for patient in sample_list.keys():
        data[patient] = {}
        for tissue in x_ticks:
            if patient not in data[patient]:
                data[patient][tissue] = []
            data[patient][tissue] = tissue_rank_dict[patient][tissue]
    plot_data = []
    group_labels = []
    subgroup_labels = []
    for group, subgroups in data.items():
        for subgroup, values in subgroups.items():
            plot_data.extend(values)
            group_labels.extend([group] * len(values))
            subgroup_labels.extend([subgroup] * len(values))
    plt.figure(figsize=(16, 7), dpi=500)
    sns.boxplot(x=subgroup_labels, y=plot_data, hue=group_labels)
    # 设置图形属性
    bwith = 3
    ax = plt.gca()  # 获取边框
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['bottom'].set_linewidth(bwith)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    plt.tick_params(axis='both', which='major', width=bwith, length=bwith * 3, labelsize=16)
    plt.tick_params(axis='x', rotation=90)
    plt.tick_params(axis='y', rotation=90)
    plt.ylabel('Ranking', size=16)
    plt.title('tissue rank')
    plt.legend(title='Dataset', bbox_to_anchor=(0, 1), loc='upper right')
    plt.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.35)
    plt.savefig(os.path.join(save_path, f'{cell} tissue in cell rank.png'))


#%%
df = pd.read_csv(single_cell_expression,header=0,index_col=0)
dict_data = df.to_dict('dict')
tissue_correction = {}
patient_type = {}
for type in sample_list:
    for patient in sample_list[type]:
        patient_type[patient] = type
#%%
for patient in patient_cell_corelation:
    type = patient_type[patient]
    tissue_correction_temp = {}
    for i in patient_cell_corelation[patient]:
        if i not in tissue_correction_temp:
            tissue_correction_temp[i] = []
        tissue_correction_temp[i].append(patient_cell_corelation[patient][i])
    tissue_correction_temp = dict(sorted(tissue_correction_temp.items(), key=lambda x: x[1],reverse=True))
    rank_dict = {}
    for index, key in enumerate(tissue_correction_temp):
        rank_dict[key] = index
    if type not in tissue_correction:
        tissue_correction[type] = {}
    for index,i in enumerate(tissue_correction_temp):
        celltype,cell, tissue = i.split('#')
        if cell not in tissue_correction[type]:
            tissue_correction[type][cell] = []
        tissue_correction[type][cell].append(index)
#%%
volcano_control = control_group_name
for volcano_compare in comapre_group_name:
    tissue_vc = {}
    for tissue in tissue_correction[volcano_control]:
        s1 = tissue_correction[volcano_control][tissue]
        s2 = tissue_correction[volcano_compare][tissue]
        fc = np.log2(np.mean(s1)/np.mean(s2))
        t_stat, p_value = stats.ttest_ind(s1, s2)
        tissue_vc[tissue] = [fc,p_value]
    # 创建示例数据

    x = []
    y_row = []
    y_log = []
    label = []
    for tissue in tissue_vc.keys():
        x.append(tissue_vc[tissue][0])
        y_row.append(tissue_vc[tissue][1])
        y_log.append(-np.log10(tissue_vc[tissue][1]))
        label.append(tissue)
    x = np.array(x)
    y_row = np.array(y_row)
    y_log = np.array(y_log)

    # 创建布尔掩码
    mask_y = y_row < 0.05
    mask_x_positive = (x > 0) & mask_y
    mask_x_negative = (x < 0) & mask_y

    # 绘制散点图
    plt.figure(figsize=(7, 7),dpi=500)
    plt.scatter(x[mask_x_positive], y_log[mask_x_positive], color='red')
    plt.scatter(x[mask_x_negative], y_log[mask_x_negative], color='blue')
    plt.scatter(x[~mask_y], y_log[~mask_y], color='gray')  # 不符合条件的点

    for i in range(len(x)):
        if mask_y[i]:
            plt.text(x[i], y_log[i], label[i], fontsize=5, ha='right' if x[i] > 0 else 'left', color='green')

    # 设置图形属性
    bwith = 3
    ax = plt.gca()  # 获取边框
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['bottom'].set_linewidth(bwith)
    plt.tick_params(axis='both', which='major', width=bwith, length=bwith*3, labelsize=16)
    plt.tick_params(axis='x', rotation=90)
    plt.xlabel(f'log2(Rank {volcano_control}./Rank {volcano_compare})',size=16)
    plt.ylabel('-log10(P-value)',size=16)
    plt.title('Cell type volcano')
    plt.tight_layout()
    plt.savefig(os.path.join(save_path,f'volcano_{volcano_control}_{volcano_compare}.png'))

#%%
volcano_control = control_group_name
for volcano_compare in comapre_group_name:
    tissue_vc = {}
    for tissue in tissue_correction[volcano_control]:
        s1 = tissue_correction[volcano_control][tissue]
        s2 = tissue_correction[volcano_compare][tissue]
        fc = np.log2(np.mean(s1)/np.mean(s2))
        t_stat, p_value = stats.ttest_ind(s1, s2)
        tissue_vc[tissue] = [fc,p_value]
    # 创建示例数据

    x = []
    y_row = []
    y_log = []
    label = []
    for tissue in tissue_vc.keys():
        x.append(tissue_vc[tissue][0])
        y_row.append(tissue_vc[tissue][1])
        y_log.append(-np.log10(tissue_vc[tissue][1]))
        label.append(tissue)
    x = np.array(x)
    y_row = np.array(y_row)
    y_log = np.array(y_log)

    # 创建布尔掩码
    mask_y = y_row < 0.05
    mask_x_positive = (x > 0) & mask_y
    mask_x_negative = (x < 0) & mask_y

    # 绘制散点图
    plt.figure(figsize=(7, 7),dpi=500)
    plt.scatter(x[mask_x_positive], y_log[mask_x_positive], color='red')
    plt.scatter(x[mask_x_negative], y_log[mask_x_negative], color='blue')
    plt.scatter(x[~mask_y], y_log[~mask_y], color='gray')  # 不符合条件的点

    # 设置图形属性
    bwith = 3
    ax = plt.gca()  # 获取边框
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['bottom'].set_linewidth(bwith)
    plt.tick_params(axis='both', which='major', width=bwith, length=bwith*3, labelsize=16)
    plt.tick_params(axis='x', rotation=90)
    plt.xlabel(f'log2(Rank {volcano_control}./Rank {volcano_compare})',size=16)
    plt.ylabel('-log10(P-value)',size=16)
    plt.title('Cell type volcano')
    plt.tight_layout()
    plt.savefig(os.path.join(save_path,f'volcano_{volcano_control}_{volcano_compare}.label.png'))