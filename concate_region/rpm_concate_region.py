#%%
import numpy as np
from sklearn.decomposition import NMF
import pandas as pd
import pysam
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import argparse
import csv
import pysam
import numpy as np
from scipy.stats import entropy
import argparse
import csv
import os

parser = argparse.ArgumentParser(description="Process samples and BED file to generate RPM results")
parser.add_argument(
    "--sample_file",
    type=str,
    required=True,
    help="Path to the sample list file (one sample per line)"
)
parser.add_argument(
    "--bed_file",
    type=str,
    required=True,
    help="Path to the BED file with genomic regions"
)
parser.add_argument(
    "--save_dir",
    type=str,
    required=True,
    help="Directory to save the output results"
)
args = parser.parse_args()

sample_file = args.sample_file
bed_file = args.bed_file
save_dir = args.save_dir


# sample_file = '/mnt/dfc_data2/project/linyusen/project/31_cfdna_wps/project/concate_region/samples.lst'
# bed_file = '/mnt/dfc_data2/project/zhoujj/project/35cfdna/04PanCancer/ref/gencode.v45.annotation.tss.b5k.bed'
# save_dir = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/concate_region/pfe'

os.makedirs(save_dir, exist_ok=True)


f = open(sample_file)
r = csv.reader(f,delimiter=' ')
sample_dict = {}
for line in r:
    group = line[0]
    name = line[1]
    path = line[2]
    if group not in sample_dict:
        sample_dict[group] = {}
    sample_dict[group][name] = path
f.close()

f = open(bed_file)
r = csv.reader(f,delimiter='\t')
region_data = []
for line in r:
    chr = line[0]
    start = int(line[1])
    end = int(line[2])
    gene = line[3]
    region_data.append([chr, start, end, gene])
f.close()
region_data = region_data[:100]
#%%

def cal_read_count(sf,chrid,start, end):
    iter = sf.fetch(chrid, start, end)
    count = 0
    for x in iter:
        count += 1
    return count


pos_dict = {}
for group in sample_dict:
    pos_dict[group] = {}
    for sample_name in sample_dict[group]:
        pos_dict[group][sample_name] = 0
        bam_file = sample_dict[group][sample_name]
        sf = pysam.AlignmentFile(bam_file, 'rb')
        for chrid, start, end, gene in tqdm(region_data,desc=sample_name):
            count = cal_read_count(sf, chrid, start, end)
            count_per_million_base = count/(abs(start-end)+1)*1000000
            pos_dict[group][sample_name] += count_per_million_base
#%%

for group in pos_dict:
    f = open(os.path.join(save_dir, group + '.rpm.tsv'),'w')
    w = csv.writer(f,delimiter='\t')
    for sample_name in pos_dict[group]:
        w.writerow([sample_name, pos_dict[group][sample_name]])
    f.close()