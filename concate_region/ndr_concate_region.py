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

parser = argparse.ArgumentParser(
        description="Extract regions from BAM file relative to TSS (center) and save results."
    )

parser.add_argument("--bam_file", required=True, help="Path to BAM file")
parser.add_argument("--bed_file", required=True, help="Path to BED file (regions of interest, e.g. TSS)")
parser.add_argument("--save_dir", required=True, help="Directory to save output files")
parser.add_argument("--sample_name", required=True, help="Sample name")
parser.add_argument("--region_length", type=int, default=2000,
                    help="Region length around center (default: 2000)")
parser.add_argument("--upstream", type=int, default=3000,
                    help="Distance upstream (5') from gene region center (default: 3000)")
parser.add_argument("--downstream", type=int, default=3000,
                    help="Distance downstream (3') from gene region center (default: 3000)")

args = parser.parse_args()

bam_file = args.bam_file
bed_file = args.bed_file
save_dir = args.save_dir
sample_name = args.sample_name
region_length = args.region_length
up_rpos = args.upstream
down_rpos = args.downstream


# bam_file = '/mnt/dfc_data2/project/zhoujj/project/ngs.data/20250731.CNGB.JWZ.cfDNA.collected/PDAC/standard/PL230613018/01alignment/PL230613018.sorted.rmdup.bam'
# bed_file = '/mnt/dfc_data2/project/zhoujj/project/35cfdna/04PanCancer/ref/gencode.v45.annotation.tss.b5k.bed'
# save_dir = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/concate_region/pfe'
# sample_name = 'PL230613018'
# region_length = 2000
# up_rpos = 3000
# down_rpos = 3000

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

def cal_coverage(sf,chrid,start, end):
    binedges = range(-1000, 1000, 5)
    pos = []
    iter = sf.fetch(chrid, start, end)  # 提取该区域的配对 reads
    tss = int(start/2 + end/2)
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
    return pos

sf = pysam.AlignmentFile(bam_file, 'rb')
med_pos = []
up_pos = []
down_pos = []
binedges = range(-int(region_length/2), int(region_length/2), 5)
for chrid, start, end, gene in tqdm(region_data):
    med_start = int((start + end) / 2) - region_length/2
    med_end = int((start + end) / 2) + region_length / 2

    up_start = int((start + end)/2 - up_rpos - region_length/2)
    up_end = int((start + end) / 2 - up_rpos + region_length/2)

    down_start = int((start + end) / 2 + down_rpos - region_length / 2)
    down_end = int((start + end) / 2 + down_rpos + region_length / 2)

    med_pos.extend(cal_coverage(sf, chrid, med_start, med_end))
    up_pos.extend(cal_coverage(sf, chrid, up_start, up_end))
    down_pos.extend(cal_coverage(sf, chrid, down_start, down_end))
#%%
med_pos = np.histogram(med_pos, binedges)  # 统计直方图（TSS ±1kb）
up_pos = np.histogram(up_pos, binedges)  # 统计直方图（TSS ±1kb）
down_pos = np.histogram(down_pos, binedges)  # 统计直方图（TSS ±1kb）
#%%
upMean = np.mean(up_pos[0])  # 平均深度
downMean = np.mean(down_pos[0])  # 平均深度
updownMean = (upMean + downMean) / 2
vectt_pos = [float(depth / updownMean) for depth in med_pos[0]]

#%%
w1 = ['sample']
w2 = [sample_name]
for index,value in enumerate(vectt_pos):
    w1.append(index)
    w2.append(value)


#%%
header_save_file = os.path.join(save_dir, 'ndr_header.tsv')
value_save_file = os.path.join(save_dir, 'ndr_value.tsv')

f = open(header_save_file,'w')
w = csv.writer(f)
w.writerow(w1)
f.close()
f = open(value_save_file,'w')
w = csv.writer(f)
w.writerow(w2)
f.close()
