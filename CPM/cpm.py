import pysam
import numpy as np
import argparse
import csv

def cpm_signal(sf,chrid,start,end,strand):
    # sf = pysam.AlignmentFile(bam_file, "rb")
    binedges = range(-1000, 1000, 1)
    s = start
    e = end
    tss = int((start+end)/2)
    tss_start = tss

    count = 0
    iter = sf.fetch(chrid, s-1000, e+1000)  # 提取该区域的配对 reads
    for x in iter:
        flaglen = abs(int(x.template_length))  # 获取插入片段长度

        if x.flag == 99 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start + flaglen
                for p in range(xs, xe):
                    rpos = p - tss  # 相对 TSS 的距离
                    if -150 < rpos < 50:
                        count+=1
                        break
        elif x.flag == 83 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start - flaglen
                for p in range(xe, xs):
                    rpos = p - tss
                    if -150 < rpos < 50:
                        count+=1
                        break


    vectt_pos = [count]
    return vectt_pos

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bed_file', type=str, required=True,
                        help='Path to the BED file defining target gene regions (e.g., upstream promoters), only support -1k, +1k regions')
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
    sf = pysam.AlignmentFile(bam_file)
    f = open(bed_file)

    lines = f.readlines()
    f.close()
    tsv_data = []
    for line in lines:
        tsv_data.append(line[:-1].split('\t'))

    cpm_data = {}
    for data_tsv in tsv_data:
        chrid, start, end, gene,strand = data_tsv[:5]
        cpm_data[gene] = cpm_signal(sf, chrid, int(start), int(end),strand)
    # 找到最长的 WPS 数组长度（用于填充）
    max_len = max(len(values) for values in cpm_data.values())
    # 保存为 CSV 文件
    with open(tsv_file, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        # 写表头：第一个是 gene 名，后面是每个位置的索引
        header = ['gene'] + [f'{i}' for i in range(max_len)]
        writer.writerow(header)
        # 写入每一行基因和对应的 WPS 值（不足部分填充为空）
        for gene, values in cpm_data.items():
            row = [gene] + list(values) + [''] * (max_len - len(values))

            writer.writerow(row)



