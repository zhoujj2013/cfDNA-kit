#%%
bam_list = {
    'A1-1': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-1/A1-1.sorted.rmdup.bam',
    'A1-2': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-2/A1-2.sorted.rmdup.bam',
    'A1-3': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-3/A1-3.sorted.rmdup.bam',
    'A1-4': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-4/A1-4.sorted.rmdup.bam',
    'A1-5-1': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-5-1/A1-5-1.sorted.rmdup.bam',
    'A1-5-2': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-5-2/A1-5-2.sorted.rmdup.bam',
    'A1-6-1': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-6-1/A1-6-1.sorted.rmdup.bam',
    'A1-6-2': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-6-2/A1-6-2.sorted.rmdup.bam',
    'A1-7': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-7/A1-7.sorted.rmdup.bam',
    'A1-8': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-8/A1-8.sorted.rmdup.bam',
    'A1-9': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-9/A1-9.sorted.rmdup.bam',
    'A1-10-1': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-10-1/A1-10-1.sorted.rmdup.bam',
    'A1-10-2': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-10-2/A1-10-2.sorted.rmdup.bam',
    'A1-11': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-11/A1-11.sorted.rmdup.bam',
    'A1-12': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/A1-12/A1-12.sorted.rmdup.bam',
    'G1-1': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-1/G1-1.sorted.rmdup.bam',
    'G1-2': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-2/G1-2.sorted.rmdup.bam',
    'G1-3': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-3/G1-3.sorted.rmdup.bam',
    'G1-4': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-4/G1-4.sorted.rmdup.bam',
    'G1-5': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-5/G1-5.sorted.rmdup.bam',
    'G1-6': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-6/G1-6.sorted.rmdup.bam',
    'G1-7': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-7/G1-7.sorted.rmdup.bam',
    'G1-8': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-8/G1-8.sorted.rmdup.bam',
    'G1-9': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-9/G1-9.sorted.rmdup.bam',
    'G1-10': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-10/G1-10.sorted.rmdup.bam',
    'G1-11': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-11/G1-11.sorted.rmdup.bam',
    'G1-12-1': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-12-1/G1-12-1.sorted.rmdup.bam',
    'G1-12-2': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-12-2/G1-12-2.sorted.rmdup.bam',
    'G1-13': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-13/G1-13.sorted.rmdup.bam',
    'G1-14': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-14/G1-14.sorted.rmdup.bam',
    'G1-15': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-15/G1-15.sorted.rmdup.bam',
    'G1-16': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-16/G1-16.sorted.rmdup.bam',
    'G1-17': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-17/G1-17.sorted.rmdup.bam',
    'G1-18': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-18/G1-18.sorted.rmdup.bam',
    'G1-19': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-19/G1-19.sorted.rmdup.bam',
    'G1-20': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-20/G1-20.sorted.rmdup.bam',
    'G1-21': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-21/G1-21.sorted.rmdup.bam',
    'G1-22': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-22/G1-22.sorted.rmdup.bam',
    'G1-7-1': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-7-1/G1-7-1.sorted.rmdup.bam',
    'G1-7-2': '/mnt/dfc_data3/project/linyusen/database/10_cfdna/01_aging/bam/G1-7-2/G1-7-2.sorted.rmdup.bam',
}
import os
cmd = '/mnt/dfc_data1/home/linyusen/miniconda/envs/lys/bin/python /mnt/dfc_data2/project/linyusen/project/31_cfdna_wps/cfDNA-kit-main/region_cpm/cpm.py --bed_file /mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/hg38.tss.1000.bed --bam_file {bam_file} --o {save_path}'
save_dir = '/mnt/dfc_data3/project/linyusen/database/10_cfdna/03_pso/region_count'
for name in bam_list:
    bam_file = bam_list[name]
    save_path = os.path.join(save_dir,name+'.csv')
    cmd1 = cmd.format(bam_file=bam_file,save_path=save_path)
    print(cmd1)