#%%
import csv
import pandas as pd
import pysam
from collections import defaultdict
import numpy as np
import itertools
import argparse


def get_five_prime_motif_from_fragment(read, n=4):
    """
    从 fragment 的 5'端提取前 n 个碱基
    read: pysam.AlignedSegment 对象（R1 或 R2）
    n: 基序长度（默认为4）
    """
    if read.is_reverse:
        # 负链：需要取反向互补
        seq = read.get_forward_sequence()[:n]
        # 取反向互补
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        motif = ''.join(complement[base] for base in seq[::-1])
    else:
        # 正链：直接取序列
        motif = read.query_sequence[:n]

    return motif
def get_fragments_in_region(bam_file, contig, start, end, n=4):
    """
    获取指定范围内的所有 fragments，并提取每个 fragment 5'端的基序

    参数:
        bam_file: BAM 文件路径
        contig: 染色体名称（如 'chr1'）
        start: 起始位置（0-based）
        end: 结束位置（0-based）
        n: 基序长度（默认4）

    返回:
        fragments: 列表，每个元素为 (fragment_name, motif, fragment_start, fragment_end, strand)
    """
    fragments = []

    # 转换为 0-based（pysam 使用 0-based）
    start = start - 1
    end = end

    # 打开 BAM 文件和参考基因组
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # 用于暂存未配对的 reads（按模板名存储）
        pending_reads = {}
        # 遍历指定范围内的所有 reads
        for read in bam.fetch(contig=contig, start=start, end=end):
            # 跳过未比对、次要比对、补充比对的 reads
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            # 只处理双端测序的 reads
            if not read.is_paired:
                print("警告: 当前 BAM 文件不是双端测序数据，无法准确获取 fragment")
                continue

            query_name = read.query_name

            # 取每个reads开头的4mer ，负链就取反
            # read1->
            #    5--------------3
            # 3--------------5
            #               <-read2

            motif = get_five_prime_motif_from_fragment(read, n)

            # 存储 fragment 信息
            fragments.append({
                'name': query_name,
                'motif': motif,
            })
    return fragments
def normalized_shannon_entropy(fragments):
    bases = ["A", "C", "G", "T"]
    k = 4
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    motif_dict = {}
    for kmer in kmers:
        motif_dict[kmer] = 0
    for data in fragments:
        if 'N' not in data['motif']:
            motif_dict[data['motif']] += 1
    counts = np.array(list(motif_dict.values()), dtype=float)

    # 全0情况
    if counts.sum() == 0:
        return 0.0

    # 转概率
    p = counts / counts.sum()

    # 去掉0避免 log2(0)
    p = p[p > 0]

    # Shannon entropy
    H = -np.sum(p * np.log2(p))

    # 最大熵
    H_max = np.log2(len(motif_dict))

    # 归一化
    H_norm = H / H_max

    return H_norm


parser = argparse.ArgumentParser()
parser.add_argument('--bed_file', type=str, required=True,
                    help='Path to the BED file defining target gene regions (e.g., upstream promoters)')
parser.add_argument('--bam_file', type=str, required=True,
                    help='Path to the cfDNA aligned BAM file')
parser.add_argument('--o', type=str, required=True,
                    help='Output path for the resulting (TSV format)')
args = parser.parse_args()

bam_file = args.bam_file
bed_file = args.bed_file
tsv_file = args.o


# bam_file = "/mnt/dfc_data3/project/linyusen/database/10_cfdna/03_pso/filter_bam/pso/PL230613012.sorted.rmdup.bam"
# bed_file = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/hg38.tss.1000.bed'
# tsv_file = '/mnt/dfc_data2/project/linyusen/project/31_cfdna_wps/cfDNA-kit-main/MDS/mds.test'
f = open(bed_file)
lines = f.readlines()
f.close()
tsv_data = []
for line in lines:
    tsv_data.append(line[:-1].split('\t'))

#%%
wps_data = {}
for line in tsv_data:
    chr = line[0]
    start = int(line[1])
    end = int(line[2])
    gene = line[3]

    fragments = get_fragments_in_region(bam_file,chr, start, end, n=4)
    entropy = normalized_shannon_entropy(fragments)
    print(gene,entropy)
    wps_data[gene] = entropy


f = open(tsv_file,'w')
w = csv.writer(f,delimiter='\t')
for gene in wps_data:
    w.writerow([gene,wps_data[gene]])
f.close()
