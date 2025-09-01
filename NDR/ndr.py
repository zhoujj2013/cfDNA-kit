import pysam
import numpy as np
import argparse
import csv

def ndr_signal(sf,chrid,start,end):
    # sf = pysam.AlignmentFile(bam_file, "rb")
    binedges = range(-1000, 1000, 5)
    s = start
    e = end
    tss = int((start+end)/2)
    tss_start = tss

    pos = []
    iter = sf.fetch(chrid, s, e)  # 提取该区域的配对 reads
    for x in iter:
        flaglen = abs(int(x.template_length))  # 获取插入片段长度

        if x.flag == 99 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start + flaglen
                for p in range(xs, xe):
                    rpos = p - tss  # 相对 TSS 的距离
                    if -1000 < rpos < 1000:
                        pos.append(rpos)
        elif x.flag == 83 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start - flaglen
                for p in range(xe, xs):
                    rpos = p - tss
                    if -1000 < rpos < 1000:
                        pos.append(rpos)

    pos.append(0)  # 避免空数组报错
    vect_pos = np.histogram(pos, binedges)  # 统计直方图（TSS ±1kb）

    ### ========= 上游（-3k ~ -1k）平均深度 =========
    utss = tss_start-2000
    us = tss_start - 3000
    ue = tss_start - 1000
    upos = []
    iter = sf.fetch(chrid, us, ue)
    for x in iter:
        flaglen = abs(int(x.template_length))
        if x.flag == 99 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start + flaglen
                for p in range(xs, xe):
                    rpos = p - utss
                    if -1000 < rpos < 1000:
                        upos.append(rpos)
        elif x.flag == 83 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start - flaglen
                for p in range(xe, xs):
                    rpos = p - utss
                    if -1000 < rpos < 1000:
                        upos.append(rpos)

    upos.append(0)
    vect_upos = np.histogram(upos, binedges)
    upMean = np.mean(vect_upos[0])  # 上游平均深度

    ### ========= 下游（+1k ~ +3k）平均深度 =========
    dtss = tss_start + 2000
    ds = tss_start + 1000
    de = tss_start + 3000
    dpos = []
    iter = sf.fetch(chrid, ds, de)
    for x in iter:
        flaglen = abs(int(x.template_length))
        if x.flag == 99 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start + flaglen
                for p in range(xs, xe):
                    rpos = p - dtss
                    if -1000 < rpos < 1000:
                        dpos.append(rpos)
        elif x.flag == 83 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start - flaglen
                for p in range(xe, xs):
                    rpos = p - dtss
                    if -1000 < rpos < 1000:
                        dpos.append(rpos)

    dpos.append(0)
    vect_dpos = np.histogram(dpos, binedges)
    downMean = np.mean(vect_dpos[0])  # 下游平均深度

    updownMean = (upMean + downMean) / 2  # 上下游平均深度
    vectt_pos = [float(depth / updownMean) for depth in vect_pos[0]]  # 归一化
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

    wps_data = {}
    for data_tsv in tsv_data:
        chrid, start, end, gene = data_tsv[:4]
        wps_data[gene] = ndr_signal(sf, chrid, int(start), int(end))
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


