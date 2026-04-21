import os
import csv

####更改的参数####
sample_file = '/mnt/dfc_data3/project/linyusen/database/10_cfdna/02_aging/sample.lst'  #sample文件
save_dir = '/mnt/dfc_data3/project/linyusen/database/10_cfdna/02_aging'  #保存目录
tmp_dir = '/mnt/dfc_data3/project/linyusen/database/10_cfdna/02_aging/temp' #缓存目录
shell_path = '/mnt/dfc_data3/project/linyusen/database/10_cfdna/02_aging/shell'  #shell目录
# scdata_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/nerous/single_cell/single_cell.expression.csv'
scdata_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/pso_hc/single_cell/single_cell.expression.csv'  #单细胞表达矩阵
bed_file = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/hg38.tss.10000.bed'  #tss位点bed文件
python = '/mnt/dfc_data1/home/linyusen/miniconda/envs/lys/bin/python' #符合运行环境的python路径
thread = 16
#####

out_dir = save_dir
if os.path.exists(out_dir) == False:
    os.makedirs(out_dir)
if os.path.exists(tmp_dir) == False:
    os.makedirs(tmp_dir)
if os.path.exists(shell_path) == False:
    os.makedirs(shell_path)
sample_list = {}
with open(sample_file,'r') as f:
    reader = csv.reader(f,delimiter='\t')
    for line in reader:
        if len(line) < 4:
            continue
        type,_,name,path = line
        print(type,_,name,path)
        if type not in sample_list:
            sample_list[type] = {}
        sample_list[type][name] = path

#%%
script = '/mnt/dfc_data2/project/linyusen/project/cfDNA_kit/pipeline/WPS/script/cfdna_analyse.py'
f_shell = open(os.path.join(shell_path,'preprocess.sh'),'w')
preprocess_path = os.path.join(out_dir,'data')
if os.path.exists(preprocess_path) == False:
    os.makedirs(preprocess_path)
for type in sample_list.keys():
    for name in sample_list[type].keys():
        print(name)
        bam_path = sample_list[type][name]
        save_path = os.path.join(preprocess_path, name)
        temp_path = os.path.join(tmp_dir, name)
        cmd = f'{python} {script} -bed_file {bed_file} -bam_file {bam_path} -out_dir {save_path} -temp_dir {temp_path} -thread {thread}\n'
        f_shell.write(cmd)
f_shell.close()
#%%
f_shell = open(os.path.join(shell_path,'summary.sh'),'w')
summary_path = os.path.join(out_dir,'summary')
if os.path.exists(summary_path) == False:
    os.makedirs(summary_path)

script = '/mnt/dfc_data2/project/linyusen/project/cfDNA_kit/pipeline/WPS/script/summary.py'
print(scdata_path)
cmd = f'{python} {script} -single_cell_matrix {scdata_path} -data_path {preprocess_path} -save_path {summary_path}\n'
f_shell.write(cmd)


script = '/mnt/dfc_data2/project/linyusen/project/cfDNA_kit/pipeline/WPS/script/figure.py'
figure_savepath = os.path.join(out_dir,'figure_198')
gene_intensity_198_path = os.path.join(summary_path,'gene_intensity_198.pkl')
patient_cell_corelation_198_path = os.path.join(summary_path,'patient_cell_corelation_198.pkl')
cmd = f'{python} {script} -gene_intensity_path {gene_intensity_198_path} -patient_cell_corelation_path {patient_cell_corelation_198_path} -single_cell_expression {scdata_path} -sample_file {sample_file} -save_path {figure_savepath}\n'
f_shell.write(cmd)
f_shell.close()