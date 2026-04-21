#%%
import numpy as np
from sklearn.decomposition import NMF
import pandas as pd
import pysam
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse

parser = argparse.ArgumentParser(description="Process BAM file and save WPS pickle")

parser.add_argument(
    '--bam_file',
    type=str,
    required=True,
    help='Path to the input BAM file'
)

parser.add_argument(
    '--save_file',
    type=str,
    required=True,
    help='Path to the output pickle file'
)

args = parser.parse_args()

bam_file = args.bam_file
save_file = args.save_file

# bam_file = '/mnt/dfc_data2/project/zhoujj/project/ngs.data/20250731.CNGB.JWZ.cfDNA.collected/leprosy/mafeng.batch2/out/DP8480005039BL_L01_251/01alignment/DP8480005039BL_L01_251.sorted.rmdup.bam'
# save_file = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/pso_nmf/data/DP8480005039BL_L01_251.pickle'

bed_file = '/mnt/dfc_data2/project/linyusen/project/31_cfdna_wps/wps_process/database/hg38.tss.10000.bed'
data_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/nc_data_single_cell/Tabula_Sapiens_All_cells_averages_assays_RNA_data.txt'

df = pd.read_csv(data_path, sep='\t')
dict_data = df.to_dict('dict')

tissue = set()
for key in dict_data:
    tissue.add(key.split('#')[-1])

tissue_single = {}
for key in dict_data:
    i = key.split('#')[-1]
    tissue_single[i] = {}
    for j in dict_data[key]:
        tissue_single[i][j] = 0
for key in dict_data:
    i = key.split('#')[-1]
    for j in dict_data[key]:
        tissue_single[i][j]+=dict_data[key][j]

tissues = sorted(tissue_single.keys())
gene_features = sorted({f for d in tissue_single.values() for f in d.keys()})


def isSoftClipped(cigar):
  for (op,count) in cigar:
    if op in [4,5,6]: return True
  return False
def aln_length(cigarlist):
  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
  return tlength
def wps_signal(bam_file,windows,chrid,start,end):
    sf = pysam.AlignmentFile(bam_file, "rb")
    '''wps'''
    length = end-start
    iter = sf.fetch(chrid, start-60-1, end+60+1)

    seq_dict = []
    for seq in iter:
        seq_dict.append(seq)
    wps_arr = np.zeros([length])
    correct = True
    def isSoftClipped(cigar):
      for (op,count) in cigar:
        if op in [4,5,6]: return True
      return False
    def aln_length(cigarlist):
      tlength = 0
      for operation,length in cigarlist:
        if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
      return tlength
    for index,read in enumerate(seq_dict):
        if read.is_duplicate or read.is_qcfail or read.is_unmapped: continue
        if isSoftClipped(read.cigar): continue
        if read.is_paired:
            if read.mate_is_unmapped: continue
            if read.rnext != read.tid: continue
            if read.is_read1 or (read.is_read2 and read.pnext + read.qlen < start - 60 - 1):
                if read.isize == 0: continue
                rstart = min(read.pos, read.pnext) + 1
                lseq = abs(read.isize)
                rend = rstart + lseq - 1
                if ((lseq < 120) or (lseq > 180)): continue

                align_start = rstart
                align_end = rend

                if align_start+windows > align_end:
                    for i in range(align_start-60,align_end+60):
                        if i - start >= 0 and i - start < length:
                            wps_arr[i - start] -= 1
                if align_start+windows < align_end:
                    for i in range(align_start+60,align_end-60):
                        if i - start >= 0 and i - start < length:
                            wps_arr[i - start] += 1
                    for i in range(align_end-60,align_end+60):
                        if i - start >= 0 and i - start < length:
                            wps_arr[i - start] -= 1
                    for i in range(align_start-60,align_start+60):
                        if i-start >= 0 and i-start < length:
                            wps_arr[i - start] -= 1
        else:
            rstart = read.pos + 1  # 1-based
            lseq = aln_length(read.cigar)
            rend = rstart + lseq - 1  # end included
            if ((lseq < 120) or (lseq > 180)): continue
            align_start = rstart
            align_end = rend

            if align_start + windows > align_end:
                for i in range(align_start - 60, align_end + 60):
                    if i - start >= 0 and i - start < length:
                        wps_arr[i - start] -= 1
            if align_start + windows < align_end:
                for i in range(align_start + 60, align_end - 60):
                    if i - start >= 0 and i - start < length:
                        wps_arr[i - start] += 1
                for i in range(align_end - 60, align_end + 60):
                    if i - start >= 0 and i - start < length:
                        wps_arr[i - start] -= 1
                for i in range(align_start - 60, align_start + 60):
                    if i - start >= 0 and i - start < length:
                        wps_arr[i - start] -= 1
    wps_arr_x = []
    for i in range(start,end):
        wps_arr_x.append(i)
    return wps_arr,wps_arr_x
def coverage_signal(bam_file,chrid,start,end):
    coverage = []
    coverage_x = []
    sf = pysam.AlignmentFile(bam_file, "rb")
    for i in range(start, end):
        coverage_x.append(i)
        coverage.append(np.sum(sf.count_coverage(chrid, i, i + 1)))
    return coverage,coverage_x
def nlms(input_signal,noisy_signal,filter_order=80,mu = 0.8,eps = 0.001 ):
    # NLMS 参数
     # 避免除零
    # 初始化
    w = np.zeros(filter_order)
    n_samples = len(input_signal)
    y = np.zeros(n_samples)
    e = np.zeros(n_samples)
    # 输入参考信号：延迟后的 noisy_signal（仿真实际中参考通道）
    x = np.concatenate((np.zeros(filter_order), input_signal))

    # NLMS 自适应滤波过程
    for n in range(n_samples):
        x_n = x[n + filter_order - 1 : n + -1 : -1]  # x[n], x[n-1], ..., x[n-L+1]
        if len(x_n) < filter_order:
            x_n = np.pad(x_n, (0, filter_order - len(x_n)))
        y[n] = np.dot(w, x_n)
        e[n] = noisy_signal[n] - y[n]
        norm = np.dot(x_n, x_n) + eps
        w += (mu / norm) * e[n] * x_n
    return y,e



rank_dict = {}

total_weight_value = {}

f = open(bed_file)
lines = f.readlines()
f.close()
tsv_data = []
for line in lines:
    tsv_data.append(line[:-1].split('\t'))

from concurrent.futures import ThreadPoolExecutor, as_completed


def process_region(chrid, start, end, gene):
    if start > end:
        start, end = end, start
    # 每个线程单独打开 BAM 文件（线程安全）
    wps_arr, wps_arr_x = wps_signal(bam_file, 120, chrid, start, end)
    return gene, wps_arr

gene_wps_dict = {}

# 多线程处理
with ThreadPoolExecutor(max_workers=4) as executor:
    futures = []
    for data in tsv_data:
        chrid = data[0]
        start = data[1]
        end = data[2]
        gene = data[3]
        start = int(start)
        end = int(end)
        start = max([start - 4000,0])
        end = end + 4000
        futures.append(executor.submit(process_region, chrid, start, end, gene))

    for future in as_completed(futures):
        result = future.result()
        if result:
            gene, wps_arr = result
            gene_wps_dict[gene] = wps_arr


import pickle
with open(save_file, "wb") as f:
    pickle.dump(gene_wps_dict, f)

