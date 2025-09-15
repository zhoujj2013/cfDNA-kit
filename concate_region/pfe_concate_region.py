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

def parse_args():
    parser = argparse.ArgumentParser(description="Process BAM file with BED regions and save results.")
    parser.add_argument(
        "--bam_file",
        type=str,
        required=True,
        help="Path to the input BAM file (e.g., PL230613018.sorted.rmdup.bam)."
    )
    parser.add_argument(
        "--bed_file",
        type=str,
        required=True,
        help="Path to the BED file containing regions of interest."
    )
    parser.add_argument(
        "--save_dir",
        type=str,
        required=True,
        help="Directory where the output results will be saved."
    )
    parser.add_argument(
        "--sample_name",
        type=str,
        required=True,
        help="Sample name used for output file naming."
    )
    parser.add_argument(
        "--windows",
        type=int,
        default=2000,
        help="Window size for region analysis (default: 2000)."
    )
    return parser.parse_args()
args = parse_args()
bam_file = args.bam_file
bed_file = args.bed_file
save_dir = args.save_dir
sample_name = args.sample_name
windows = args.windows

# bam_file = '/mnt/dfc_data2/project/zhoujj/project/ngs.data/20250731.CNGB.JWZ.cfDNA.collected/PDAC/standard/PL230613018/01alignment/PL230613018.sorted.rmdup.bam'
# bed_file = '/mnt/dfc_data2/project/zhoujj/project/35cfdna/04PanCancer/ref/gencode.v45.annotation.tss.b5k.bed'
# save_dir = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/concate_region/pfe'
# sample_name = 'PL230613018'
# windows = 2000

os.makedirs(save_dir, exist_ok=True)


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
sf = pysam.AlignmentFile(bam_file)
for chrid, start, end, gene in tqdm(region_data):
    regions = split_region(start, end, windows)
    bins = list(range(50, 400, 50))

    for index,(sub_start, sub_end) in enumerate(regions):
        if index not in lengths_dict:
            lengths_dict[index] = []
        lengths_dict[index].extend(extract_fragment_lengths(sf, chrid, sub_start, sub_end))
entropy_dict = {}
for index in lengths_dict:
    bins = list(range(50, 400, 50))
    entropy_dict[index] = compute_entropy(lengths_dict[index],bins)
#%%
w1 = ['sample']
for index in entropy_dict:
    w1.append(index)
w2 = [sample_name]
for index in entropy_dict:
    w2.append(index)

header_save_file = os.path.join(save_dir, 'pfe_header.tsv')
value_save_file = os.path.join(save_dir, 'pfe_value.tsv')
f = open(header_save_file,'w')
w = csv.writer(f)
w.writerow(w1)
f.close()
f = open(value_save_file,'w')
w = csv.writer(f)
w.writerow(w2)
f.close()
