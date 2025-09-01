import pysam
import numpy as np
from scipy.stats import entropy
import argparse
import csv

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
    if len(filtered_lengths) < 20:
        return None  # too few fragments
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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bed_file', type=str, required=True,
                        help='Path to the BED file defining target gene regions (e.g., upstream promoters)')
    parser.add_argument('--bam_file', type=str, required=True,
                        help='Path to the cfDNA aligned BAM file')
    parser.add_argument('--o', type=str, required=True,
                        help='Output path for the resulting (TSV format)')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    bed_file = args.bed_file
    bam_file = args.bam_file
    tsv_file = args.o
    windows = 1000
    sf = pysam.AlignmentFile(bam_file)
    f = open(bed_file)

    lines = f.readlines()
    f.close()
    tsv_data = []
    for line in lines:
        tsv_data.append(line[:-1].split('\t'))

    wps_data = {}
    for chrid, start, end, gene in tsv_data[:4]:
        wps_data[gene] = pfe_signal(sf, windows, chrid, int(start), int(end))
    # 找到最长的 WPS 数组长度（用于填充）
    max_len = max(len(values) for values in wps_data.values())
    # 保存为 CSV 文件
    with open(tsv_file, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        # 写表头：第一个是 gene 名，后面是每个位置的索引
        header = ['gene'] + [f'{i}' for i in range(max_len)]
        writer.writerow(header)
        # 写入每一行基因和对应的 WPS 值（不足部分填充为空）
        for gene, values in wps_data.items():
            row = [gene] + list(values) + [''] * (max_len - len(values))

            writer.writerow(row)
