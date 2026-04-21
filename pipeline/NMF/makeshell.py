#%%
import os


#疾病组
aging_bam_list = [
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-1/G1-1.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-2/G1-2.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-3/G1-3.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-4/G1-4.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-5/G1-5.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-6/G1-6.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-7/G1-7.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-8/G1-8.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-9/G1-9.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-10/G1-10.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-11/G1-11.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-12-1/G1-12-1.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-12-2/G1-12-2.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-13/G1-13.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-14/G1-14.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-15/G1-15.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-16/G1-16.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-17/G1-17.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-18/G1-18.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-19/G1-19.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-20/G1-20.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-21/G1-21.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-22/G1-22.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-7-1/G1-7-1.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-7-2/G1-7-2.sorted.rmdup.bam',
]

#健康对照组
young_bam_list = [
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/PL230613001/PL230613001.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/PL230613002/PL230613002.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/PL230613003/PL230613003.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/PL230613004/PL230613004.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/PL230613005/PL230613005.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/PL230613006/PL230613006.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/PL230613007/PL230613007.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/PL230613008/PL230613008.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/PL230613009/PL230613009.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/PL230613010/PL230613010.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/PL230613012/PL230613012.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/PL230613011/PL230613011.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-1/A1-1.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-2/A1-2.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-3/A1-3.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-4/A1-4.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-5-1/A1-5-1.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-5-2/A1-5-2.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-6-1/A1-6-1.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-6-2/A1-6-2.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-7/A1-7.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-8/A1-8.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-9/A1-9.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-10-1/A1-10-1.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-10-2/A1-10-2.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-11/A1-11.sorted.rmdup.bam',
'/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-12/A1-12.sorted.rmdup.bam',
]


save_path = '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/nmf/'

shell_dir = os.path.join('/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/nmf/','shell')
os.makedirs(shell_dir,exist_ok=True)

f = open(os.path.join(shell_dir,'work1_cal_wps.sh'),'w')
wps_save_dir = os.path.join(save_path,'wps')
os.makedirs(wps_save_dir,exist_ok=True)
wps_pickle_bam = {}
for bam_file in aging_bam_list:
    bam_name = os.path.basename(bam_file)
    pl = bam_file.split('/')[-1].split('.')[0]
    save_dir = os.path.join(wps_save_dir,pl+'.pickle')
    wps_pickle_bam[pl+'.pickle'] = bam_name
    cmd = f'/mnt/dfc_data1/home/linyusen/miniconda/envs/lys/bin/python /mnt/dfc_data2/project/linyusen/project/cfDNA_kit/pipeline/NMF/script/cal_wps.py --bam_file {bam_file} --save_file {save_dir}'
    ccmd = f'qsub -o ./ -e ./ -b y -cwd  "{cmd}"\n'
    f.write(ccmd)
for bam_file in young_bam_list:
    bam_name = os.path.basename(bam_file)
    pl = bam_file.split('/')[-1].split('.')[0]
    save_dir = os.path.join(wps_save_dir,pl+'.pickle')
    wps_pickle_bam[pl+'.pickle'] = bam_name
    cmd = f'/mnt/dfc_data1/home/linyusen/miniconda/envs/lys/bin/python /mnt/dfc_data2/project/linyusen/project/cfDNA_kit/pipeline/NMF/script/cal_wps.py --bam_file {bam_file} --save_file {save_dir}'
    ccmd = f'qsub -o ./ -e ./ -b y -cwd   "{cmd}"\n'
    f.write(ccmd)
f.close()
print(f'首先在头结点运行 {os.path.join(shell_dir,"work1_cal_wps.sh")}')

import pysam
def count_reads_idxstats(bam_file):
    """
    返回 BAM 文件的 reads 总数（mapped + unmapped）
    """
    stats = pysam.idxstats(bam_file)
    total_reads = 0
    for line in stats.strip().split('\n'):
        fields = line.split('\t')
        mapped = int(fields[2])
        unmapped = int(fields[3])
        total_reads += mapped + unmapped
    return total_reads
bam_scale = {}
for bam_path in aging_bam_list:
    bam_name = os.path.basename(bam_path)
    total = count_reads_idxstats(bam_path)
    bam_scale[bam_name] = total/20000000
for bam_path in young_bam_list:
    bam_name = os.path.basename(bam_path)
    total = count_reads_idxstats(bam_path)
    bam_scale[bam_name] = total/20000000
f = open(os.path.join(shell_dir,'work2_cal_nmf.sh'),'w')
nmf_save_dir = os.path.join(save_path,'nmf')
os.makedirs(nmf_save_dir,exist_ok=True)

patirnt_weight_path_dict = {}
for i in wps_pickle_bam:
    bam_name = wps_pickle_bam[i]
    scale = bam_scale[bam_name]
    cmd = f"/mnt/dfc_data1/home/linyusen/miniconda/envs/lys/bin/python /mnt/dfc_data2/project/linyusen/project/cfDNA_kit/pipeline/NMF/script/cal_nmf.py --pickle_dir {os.path.join(wps_save_dir,i)} --csv_file {os.path.join(nmf_save_dir,i.split('.')[0]+'.weight.csv')} --scale_score {scale}"
    ccmd = f'qsub -o ./ -e ./ -b y -cwd "{cmd}"\n'
    f.write(cmd+'\n')

    for bam_path in aging_bam_list:
        if bam_name in bam_path:
            patirnt_weight_path = os.path.join(nmf_save_dir, i.split('.')[0] + '.weight.csv')
            patient_name = bam_name.split('.')[0]
            patirnt_weight_path_dict[patient_name] = patirnt_weight_path
            break
f.close()
print(f'上一步彻底运行完后运行 {os.path.join(shell_dir,"work2_cal_nmf.sh")}')

f = open(os.path.join(nmf_save_dir,'compare.lst'),'w')
compare_lst_path = os.path.join(nmf_save_dir,'compare.lst')
for bam_path in young_bam_list:
    bam_name = os.path.basename(bam_path)
    weight_path = os.path.join(nmf_save_dir, bam_name.split('.')[0] + '.weight.csv')
    f.write(weight_path+'\n')
f.close()


import os
f = open(os.path.join(shell_dir,'work3_cal_contribution.sh'),'w')
contribution_dir = os.path.join(save_path,'contribution')
os.makedirs(contribution_dir,exist_ok=True)
for patient in patirnt_weight_path_dict:
    weight_path = patirnt_weight_path_dict[patient]
    contribution = os.path.join(contribution_dir,f'{patient}.normalization.contribution.csv')
    cmd = f"/mnt/dfc_data1/home/linyusen/miniconda/envs/lys/bin/python /mnt/dfc_data2/project/linyusen/project/cfDNA_kit/pipeline/NMF/script/contribution_matrix_our.py --compare_weight_lst {compare_lst_path} --patient_weight {weight_path} --save_csv {contribution}"
    ccmd = f'qsub -o ./ -e ./ -v -b y -cwd "{cmd}"\n'
    f.write(cmd+'\n')
f.close()
print(f'上一步彻底运行完后运行 {os.path.join(shell_dir,"work3_cal_contribution.sh")}')