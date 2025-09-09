import sys
import pysam
from tqdm import tqdm
import numpy as np
from scipy import signal
import pywt
import pickle
import pandas as pd
import argparse
import csv


def parse_args():
    parser = argparse.ArgumentParser(description="""
cfDNA Signal Energy Analysis Module:
This program enhances the fragment distribution signal of cell-free DNA (cfDNA) 
within specified genomic regions. It applies adaptive filtering to highlight 
significant differential features, followed by a wavelet transform to compute 
energy in the Level-3 frequency band. The extracted energy signals can aid in 
the downstream analysis of gene expression status.
""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--gene_region_file', type=str, required=True,
                        help='Path to the gene region BED file')
    parser.add_argument('--bam_file', type=str, required=True,
                        help='Path to the cfDNA BAM file')
    parser.add_argument('--save_path', type=str, required=True,
                        help='Path to save the output energy CSV file')
    parser.add_argument('--ref_chrid', type=str, default='chr7',
                        help='Reference chromosome ID (e.g., chr7)')
    parser.add_argument('--ref_start', type=int, default=73006144,
                        help='Start position of the genomic region (e.g., 73006144)')
    parser.add_argument('--ref_end', type=int, default=73016144,
                        help='End position of the genomic region (e.g., 73016144)')
    parser.add_argument('--windows', type=int, default=120,
                        help='Sliding window size for energy computation (default: 120 bp)')

    return parser.parse_args()


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



args = parse_args()
gene_region_file = args.gene_region_file
bam_file = args.bam_file
save_path = args.save_path
ref_chrid = args.ref_chrid
ref_start = args.ref_start
ref_end = args.ref_end
windows = args.windows

# gene_region_file = '/mnt/dfc_data2/project/linyusen/database/02_hg/hg38/gene.info.bed'
# bam_file = '/mnt/dfc_data2/project/zhoujj/project/35cfdna/00baseline/cohort2_out/PL230613002/01alignment/PL230613002.sorted.rmdup.bam'
# save_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/nmf_data/PL230613002.energies.csv'
# ref_chrid = 'chr7'
# ref_start = 73006144
# ref_end = 73016144
# windows = 120

# bam_file = sys.argv[1]
# save_path = sys.argv[2]


gene_info = {}
f = open(gene_region_file)
r = csv.reader(f,delimiter='\t')
for line in r:
    chrid = line[0]
    start = int(line[1])
    end = int(line[2])
    gene = line[3]
    dir = line[4]
    gene_info[gene] = [chrid,start,end]


fs = 1000
cutoff_frequency = 100
order = 10
b, a = signal.butter(order, cutoff_frequency, fs=fs, btype='low')



energies_dict = {}
for gene in tqdm(gene_info):
    chrid, start, end = gene_info[gene]
    temp = (ref_end - ref_start) - (end - start)
    ref_start_temp = ref_start - temp
    ref_end_temp = ref_end
    ref_wps_arr, ref_wps_arr_x = wps_signal(bam_file, windows, ref_chrid, ref_start_temp, ref_end_temp)
    ref_wps_arr = signal.filtfilt(b, a, ref_wps_arr)
    signal_arr, wps_arr_x = wps_signal(bam_file, 120, chrid, start, end)
    signal_arr = signal.filtfilt(b, a, signal_arr)
    mini_length = min([len(ref_wps_arr), len(signal_arr)])
    ref_wps_arr = ref_wps_arr[:mini_length]
    signal_arr = signal_arr[:mini_length]
    filter_signal, filter_error = nlms(signal_arr, ref_wps_arr, filter_order=80, mu=0.01, eps=0.1) # Use Normalized Least Mean Squares (NLMS) adaptive filter to extract differential signal
    coeffs = pywt.wavedec(filter_error, 'db4', level=3) # Perform 3-level discrete wavelet transform (DWT) on the residual error signal
    energies = [np.sum(np.square(c)) for c in coeffs] # Calculate energy at each decomposition level (sum of squared coefficients)
    energies_dict[gene] = {}
    for index in range(len(energies)):
        energies_dict[gene][f"db4_level{index}"] = energies[index]

energies_pd = pd.DataFrame(energies_dict).T
energies_pd.to_csv(save_path)

