
import pysam
import numpy as np
import argparse
import csv
from scipy import signal

def isSoftClipped(cigar):
  for (op,count) in cigar:
    if op in [4,5,6]: return True
  return False
def aln_length(cigarlist):
  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
  return tlength
def wps_signal(sf,windows,chrid,start,end):
    # sf = pysam.AlignmentFile(bam_file, "rb")
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
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bed_file', type=str, required=True,
                        help='Path to the BED file defining target gene regions (e.g., upstream promoters)')
    parser.add_argument('--bam_file', type=str, required=True,
                        help='Path to the cfDNA aligned BAM file')
    parser.add_argument('--o', type=str, required=True,
                        help='Output path for the resulting (TSV format)')
    parser.add_argument('--ref_chrid', type=str, default='chr7',
                        help='Reference chromosome ID (e.g., chr7)')
    parser.add_argument('--ref_start', type=int, default=73006144,
                        help='Start position of the genomic region (e.g., 73006144)')
    parser.add_argument('--ref_end', type=int, default=73016144,
                        help='End position of the genomic region (e.g., 73016144)')

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    bed_file = args.bed_file
    bam_file = args.bam_file
    tsv_file = args.o
    ref_chrid = args.ref_chrid
    ref_start = args.ref_start
    ref_end = args.ref_end



    windows = 120
    fs = 1000
    cutoff_frequency = 100
    order = 10
    b, a = signal.butter(order, cutoff_frequency, fs=fs, btype='low')

    sf = pysam.AlignmentFile(bam_file)

    ref_wps_arr, ref_wps_arr_x = wps_signal(sf, windows, ref_chrid, ref_start, ref_end)
    ref_wps_arr = signal.filtfilt(b, a, ref_wps_arr)

    f = open(bed_file)

    lines = f.readlines()
    f.close()
    tsv_data = []
    for line in lines:
        tsv_data.append(line[:-1].split('\t'))

    anps_data = {}
    for chrid, start, end, gene, _ in tsv_data[:]:

        signal_arr, wps_arr_x = wps_signal(sf,windows,chrid,int(start),int(end))
        signal_arr = signal.filtfilt(b, a, signal_arr)
        filter_signal, filter_error = nlms(signal_arr, ref_wps_arr, filter_order=80, mu=0.01,eps=0.1)
        anps_data[gene] = filter_error

    # 找到最长的 WPS 数组长度（用于填充）
    max_len = max(len(values) for values in anps_data.values())
    # 保存为 CSV 文件
    with open(tsv_file, mode='w', newline='') as file:
        writer = csv.writer(file,delimiter='\t')
        # 写表头：第一个是 gene 名，后面是每个位置的索引
        header = ['gene'] + [f'{i}' for i in range(max_len)]
        writer.writerow(header)
        # 写入每一行基因和对应的 WPS 值（不足部分填充为空）
        for gene, values in anps_data.items():
            row = [gene] + list(values) + [''] * (max_len - len(values))
            writer.writerow(row)