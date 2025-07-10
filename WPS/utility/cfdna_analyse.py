import multiprocessing
import os
import argparse
import shutil
import csv
import yaml

parse = argparse.ArgumentParser()
parse.add_argument('-bed_file',required=True,help='需要计算的region bed的参考文件/mnt/dfc_data2/project/linyusen/database/46_cfdna/cfDNA-master/expression/transcriptAnno-hg38.75.upstream.bed')
parse.add_argument('-bam_file',required=True,help='cfdna的bam文件')
parse.add_argument('-out_dir',required=True,help='输出目录')
parse.add_argument('-temp_dir',required=True,help='临时文件夹，存放中间数据')
parse.add_argument('-thread',type=int,default=32,help='多线程')
parse.add_argument('-config',type=str,help='config file')
parse.add_argument('-minInsert',type=int,default=120,help='最小fragment 长度')
parse.add_argument('-maxInsert',type=int,default=180,help='最大fragment 长度')
args = parse.parse_args()
bed_file = args.bed_file
bam_file = args.bam_file
out_dir = args.out_dir
temp_dir = args.temp_dir
thread = args.thread
config = args.config
minInsert=args.minInsert
maxInsert=args.maxInsert

#%%
with open(file=config, mode="r",encoding='utf-8') as f:
    res = yaml.safe_load(f)
fft_script = res['fft_script']
wps_script = res['wps_script']
rscript_software = res['rscript_software']
python27_software = res['python27_software']
python38_software = res['python38_software']
concatfft = res['concatfft']
concatwps = res['concatwps']


if os.path.exists(temp_dir) == False:
    os.makedirs(temp_dir)
f = open(bed_file)
lines = f.readlines()
f.close()

tsv_data = []
for line in lines:
    tsv_data.append(line[:-1].split('\t'))

def split_list(lst, n):
    """
    将列表 lst 平均分成 n 份
    """
    k, m = divmod(len(lst), n)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]

def cmd_process(cmd):
    os.system(cmd)

def multiprocessing_cmd(cmd_list,tips):
    processes = []
    for cmd in cmd_list:
        p = multiprocessing.Process(target=cmd_process, args=(cmd,))
        processes.append(p)
        p.start()
    for p in processes:
        p.join()
    print('{} completed'.format(tips))
splitted_tsv_data = split_list(tsv_data, thread)
if os.path.exists(os.path.join(temp_dir,'tsv')) == False:
    os.makedirs(os.path.join(temp_dir,'tsv'))
for index,tsv in enumerate(splitted_tsv_data):
    f = open(os.path.join(temp_dir,'tsv','{}.tsv'.format(index)),'w')
    writer = csv.writer(f,delimiter='\t')
    for chr,start,end,cid,l in tsv:
        chr = chr.replace('chr','')
        writer.writerow([cid,chr,start,end,l])
    f.close()
#%%

print('wps calcualting')
if os.path.exists(os.path.join(out_dir,'wps')) == False:
    os.makedirs(os.path.join(out_dir,'wps'))
wps_out_dir = "'{}/block_%s.tsv.gz'".format(os.path.join(out_dir,'wps'))
wps_cmd_list = []

# # 使用 pysam.AlignmentFile 打开 BAM 文件
# with pysam.AlignmentFile(bam_file, "rb") as bam:
#     # 统计 BAM 文件中 reads 的数量
#     read_count = bam.count()
#     print('read_count:',read_count)

for file in os.listdir(os.path.join(temp_dir,'tsv')):
    region_file = os.path.join(temp_dir,'tsv',file)
    wps_cmd= f"{python27_software} {wps_script} --minInsert={minInsert} --maxInsert={maxInsert} -i {region_file} -o {wps_out_dir} {bam_file}"
    print(wps_cmd)
    wps_cmd_list.append(wps_cmd)
multiprocessing_cmd(wps_cmd_list,'wps')
#%%
print('fft calcualting')
if os.path.exists(os.path.join(out_dir,'fft')) == False:
    os.makedirs(os.path.join(out_dir,'fft'))
fft_out_dir = os.path.join(out_dir,'fft')
wps_out_dir = os.path.join(out_dir,'wps')
file_list = os.listdir(wps_out_dir)
splitted_file_list = split_list(file_list,thread)
if os.path.exists(os.path.join(temp_dir,'fft')) == False:
    os.makedirs(os.path.join(temp_dir,'fft'))
fft_sh_cmd= []
for index,file_list in enumerate(splitted_file_list):
    f = open(os.path.join(temp_dir, 'fft', '{}.sh'.format(index)), 'w')
    for file in file_list:
        fft_cmd = f'{rscript_software} {fft_script} {wps_out_dir} {fft_out_dir} {file}'
        f.writelines(fft_cmd+'\n')
    f.close()
    fft_sh_cmd.append('/usr/bin/bash '+os.path.join(temp_dir, 'fft', '{}.sh'.format(index)))
multiprocessing_cmd(fft_sh_cmd,'fft')

#%%

fft_save = os.path.join(out_dir,'fft.pickle')
print(fft_save)
cmd = f'{python38_software} {concatfft} -i {fft_out_dir} -s {fft_save}'
print(cmd)
os.system(cmd)


wps_save = os.path.join(out_dir,'wps.pickle')
print(wps_save)
cmd = f'{python38_software} {concatwps} -i {wps_out_dir} -s {wps_save}'
print(cmd)
os.system(cmd)

shutil.rmtree(fft_out_dir)
shutil.rmtree(wps_out_dir)