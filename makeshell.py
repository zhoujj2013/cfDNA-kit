import os
import csv
import yaml



sample_file = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/plasma_100_test/sample.lst'
save_dir = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/plasma_100_test'
tmp_dir = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/plasma_100_test/temp'
shell_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/plasma_100_test/shell'
config = '/mnt/dfc_data2/project/linyusen/project/cfdna_r_github/config/config.yaml'
thread = 70


with open(file=config, mode="r",encoding='utf-8') as f:
    res = yaml.safe_load(f)

#%%
scdata_path = res['scdata_path']
bed_file = res['bed_file']
python38_software = res['python38_software']
cfdna_analyse_script = res['cfdna_analyse_script']
summary_script = res['summary_script']
figure_script = res['figure_script']



#%%
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
        if type not in sample_list:
            sample_list[type] = {}
        sample_list[type][name] = path
#%%

f_shell = open(os.path.join(shell_path,'preprocess.sh'),'w')
preprocess_path = os.path.join(out_dir,'data')
if os.path.exists(preprocess_path) == False:
    os.makedirs(preprocess_path)
for type in sample_list.keys():
    for name in sample_list[type].keys():
        bam_path = sample_list[type][name]
        save_path = os.path.join(preprocess_path, name)
        temp_path = os.path.join(tmp_dir, name)
        cmd = f'{python38_software} {cfdna_analyse_script} -bed_file {bed_file} -bam_file {bam_path} -out_dir {save_path} -temp_dir {temp_path} -thread {thread} -config {config}\n'
        f_shell.write(cmd)
f_shell.close()
#%%
f_shell = open(os.path.join(shell_path,'summary.sh'),'w')
summary_path = os.path.join(out_dir,'summary')
if os.path.exists(summary_path) == False:
    os.makedirs(summary_path)


cmd = f'{python38_software} {summary_script} -single_cell_matrix {scdata_path} -data_path {preprocess_path} -save_path {summary_path}\n'
f_shell.write(cmd)


figure_savepath = os.path.join(out_dir,'figure_ratio')
gene_intensity_ratio_path = os.path.join(summary_path,'gene_intensity_ratio.pkl')
patient_cell_corelation_ratio_path = os.path.join(summary_path,'patient_cell_corelation_ratio.pkl')
cmd = f'{python38_software} {figure_script} -gene_intensity_path {gene_intensity_ratio_path} -patient_cell_corelation_path {patient_cell_corelation_ratio_path} -single_cell_expression {scdata_path} -sample_file {sample_file} -save_path {figure_savepath}\n'
f_shell.write(cmd)


figure_savepath = os.path.join(out_dir,'figure_198')
gene_intensity_198_path = os.path.join(summary_path,'gene_intensity_198.pkl')
patient_cell_corelation_198_path = os.path.join(summary_path,'patient_cell_corelation_198.pkl')
cmd = f'{python38_software} {figure_script} -gene_intensity_path {gene_intensity_198_path} -patient_cell_corelation_path {patient_cell_corelation_198_path} -single_cell_expression {scdata_path} -sample_file {sample_file} -save_path {figure_savepath}\n'
f_shell.write(cmd)
f_shell.close()