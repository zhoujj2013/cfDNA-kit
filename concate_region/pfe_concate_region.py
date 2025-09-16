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

parser = argparse.ArgumentParser(description="Process samples and BED file to generate PFE results")
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
def extract_fragment_lengths(bamfile, chrom, start, end):
    """Extract fragment lengths from properly paired reads within region."""
    lengths = []
    for read in bamfile.fetch(chrom, start, end):
        if read.flag in (99, 83) and read.reference_id == read.next_reference_id:
            if abs(read.template_length) < 2000:
                lengths.append(abs(read.template_length))
    return lengths
def split_region(start, end, window_size):
    """
    将给定区域 [start, end) 按 window_size 划分多个子区域。

    返回：[(start1, end1), (start2, end2), ...]
    """
    regions = []
    for i in range(start, end, window_size):
        sub_start = i
        sub_end = min(i + window_size, end)
        regions.append([sub_start, sub_end])
    return regions
def compute_entropy(lengths, bins):
    """Compute histogram and entropy of fragment length distribution."""
    filtered_lengths = [l for l in lengths if 50 <= l <= 400]
    hist, _ = np.histogram(filtered_lengths, bins=bins)
    probs = hist / np.sum(hist)
    return entropy(probs)
def pfe_signal(sf, windows, chrid, start, end):
    # bamfile = pysam.AlignmentFile(bam_file, "rb")
    regions = split_region(start,end,windows)
    bins = list(range(50, 400, 50))
    entropy_vals = []
    for sub_start, sub_end in regions:
        lengths = extract_fragment_lengths(sf, chrid, sub_start, sub_end)
        entropy_val = compute_entropy(lengths,bins)
        entropy_vals.append(entropy_val)
    return entropy_vals

lengths_dict = {}
for group in sample_dict:
    lengths_dict[group] = {}
    for sample_name in sample_dict[group]:
        lengths_dict[group][sample_name] = []
        bam_file = sample_dict[group][sample_name]
        sf = pysam.AlignmentFile(bam_file)
        for chrid, start, end, gene in tqdm(region_data,desc=sample_name):
            lengths_dict[group][sample_name].extend(extract_fragment_lengths(sf, chrid, start, end))
#%%
bins = list(range(50, 400, 50))
entropy_dict = {}
for group in sample_dict:
    entropy_dict[group] = {}
    for sample_name in sample_dict[group]:
        entropy_dict[group][sample_name] = compute_entropy(lengths_dict[group][sample_name],bins)
#%%
for group in entropy_dict:
    f = open(os.path.join(save_dir, group + '.pfe.tsv'),'w')
    w = csv.writer(f,delimiter='\t')
    for sample_name in entropy_dict[group]:
        w.writerow([sample_name, entropy_dict[group][sample_name]])
    f.close()