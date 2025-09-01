
import pysam
import numpy as np
import argparse
import csv

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
    windows = 120
    sf = pysam.AlignmentFile(bam_file)
    f = open(bed_file)

    lines = f.readlines()
    f.close()
    tsv_data = []
    for line in lines:
        tsv_data.append(line[:-1].split('\t'))

    wps_data = {}
    for data_tsv in tsv_data:
        chrid, start, end, gene = data_tsv[:4]
        wps_data[gene],_ = wps_signal(sf,windows,chrid,int(start),int(end))
        print(wps_data[gene])

    # 找到最长的 WPS 数组长度（用于填充）
    max_len = max(len(values) for values in wps_data.values())
    # 保存为 CSV 文件
    with open(tsv_file, mode='w', newline='') as file:
        writer = csv.writer(file,delimiter='\t')
        # 写表头：第一个是 gene 名，后面是每个位置的索引
        header = ['gene'] + [f'{i}' for i in range(max_len)]
        writer.writerow(header)
        # 写入每一行基因和对应的 WPS 值（不足部分填充为空）
        for gene, values in wps_data.items():
            row = [gene] + list(values) + [''] * (max_len - len(values))

            writer.writerow(row)

