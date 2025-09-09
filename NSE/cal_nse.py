#%%
import os
import pickle
import numpy as np
from scipy import signal
from scipy.signal import find_peaks
from scipy.stats import entropy
import argparse
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
from scipy.signal import find_peaks
from scipy import signal
from scipy.stats import entropy
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.spatial.distance import jensenshannon
import pysam
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Process BAM and BED files for cfDNA analysis.")

    parser.add_argument('--bed_file', type=str, required=True,
                        help="Path to the BED file, e.g., gencode.v45.annotation.tss.b5k.bed")
    parser.add_argument('--bam_file', type=str, required=True,
                        help="Path to the BAM file, e.g., PL230613013.sorted.rmdup.bam")
    parser.add_argument('--save_path', type=str, required=True,
                        help="Directory to save output files")

    args = parser.parse_args()
    return args


args = parse_args()
bed_file = args.bed_file

bam_file = args.bam_file
save_path = args.save_path

os.makedirs(save_path, exist_ok=True)
# bed_file = '/mnt/dfc_data2/project/zhoujj/project/35cfdna/04PanCancer/ref/gencode.v45.annotation.tss.b5k.bed'
# bam_file = '/mnt/dfc_data2/project/zhoujj/project/ngs.data/20250731.CNGB.JWZ.cfDNA.collected/PDAC/standard/PL230613013/01alignment/PL230613013.sorted.rmdup.bam'
# save_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/'

gene_dict = {}
f = open(bed_file)
for line in f.readlines():
    line = line.split('\t')
    chrid = line[0]
    start = min([int(line[1]),int(line[2])])
    end = max([int(line[1]),int(line[2])])
    gene = line[3]
    gene_dict[gene]=[chrid,start,end,gene]
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
def js_distance_from_samples(samples1, samples2, bins="auto"):
    # 计算 Jensen-Shannon 距离
    return jensenshannon(samples1, samples2)
def calculate_entropy(data,reference_data, bins=50):
    """
    计算熵：先对数据进行分箱，再计算每个箱的概率分布，然后用scipy计算熵
    """
    if len(data) == 0:
        return None  # 如果数据为空，返回None
    hist, _ = np.histogram(data, bins=bins, density=True)  # 直方图概率分布
    hist = hist[hist > 0]  # 去掉0避免log问题

    reference_hist, _ = np.histogram(reference_data, bins=bins, density=True)  # 直方图概率分布
    reference_hist = reference_hist[reference_hist > 0]  # 去掉0避免log问题
    return entropy(pk=hist,qk=reference_hist)

from tqdm import tqdm
signal_dict = {}
for gene in tqdm(gene_dict):
    chrid, start, end,_ = gene_dict[gene]
    wps_arr, wps_arr_x = wps_signal(bam_file, 120, chrid, start, end)
    signal_dict[gene] = wps_arr
#%%

fs = 1000  # 采样率
cutoff_frequency = 8  # 截止频率
order = 5  # 滤波器阶数
b, a = signal.butter(order, cutoff_frequency, fs=fs, btype='low')

js_distance_high = []
f1 = open(os.path.join(save_path,'hk_gene.js_distance.tsv'),'w')
w1 = csv.writer(f1,delimiter='\t')
w1.writerow(['gene','js_distance_nucleu_loc' , 'js_distance_peak_distance_value' , 'js_distance_peak_bottom_distance_value' , 'js_slope_up' , 'js_slope_down'])
f2 = open(os.path.join(save_path,'hk_gene.peak.tsv'),'w')
w2 = csv.writer(f2,delimiter='\t')
w2.writerow(['gene','q1','q3'])
f3 = open(os.path.join(save_path,'hk_gene.bottom.tsv'),'w')
w3 = csv.writer(f3,delimiter='\t')
w3.writerow(['gene','q1','q3'])
f4 = open(os.path.join(save_path,'hk_gene.nuclusome.tsv'),'w')
w4 = csv.writer(f4,delimiter='\t')
f5 = open(os.path.join(save_path,'hk_gene.peak-bottom.tsv'),'w')
w5 = csv.writer(f5,delimiter='\t')
f6 = open(os.path.join(save_path,'hk_gene.nuclusome_distance.tsv'),'w')
w6 = csv.writer(f6,delimiter='\t')
for gene in tqdm(signal_dict):
    try:
        nucleu_loc = []
        peak_value = []
        bottom_value = []
        peak_distance = []
        peak_bottom_value = []
        peak_loc = []
        bottom_loc = []
        k1_list = []
        k2_list = []
        reference_peak_value = []
        reference_bottom_value = []
        reference_nucleu_loc = []
        reference_peak_distance = []
        reference_peak_bottom_value = []
        reference_k1_list = []
        reference_k2_list = []

        wps_arr = []
        for pos in signal_dict[gene]:
            wps_arr = signal_dict[gene]
        wps_arr = np.array(wps_arr)
        filtered_signal = signal.filtfilt(b, a, wps_arr)
        filtered_signal = filtered_signal - np.mean(filtered_signal)
        q1 = np.percentile(filtered_signal, 25)  # 25% 分位
        q3 = np.percentile(filtered_signal, 75)  # 75% 分位
        w_js_distance = []
        w_nucleu_loc = [gene]
        w_peak_value = [gene, q1, q3]
        w_bottom_value = [gene, q1, q3]
        w_peak_bottom_value = [gene]
        w_peak_distance = [gene]
        peaks, _ = find_peaks(filtered_signal)
        inverted_peaks, _ = find_peaks(-filtered_signal)
        nucleusome_loc = []
        for peak in peaks:
            minima_before = inverted_peaks[inverted_peaks < peak]
            minima_after = inverted_peaks[inverted_peaks > peak]
            if len(minima_before) > 0 and len(minima_after) > 0:
                min_before = minima_before[-1]
                min_after = minima_after[0]
                if abs(filtered_signal[min_after] - filtered_signal[peak]) < 3 or abs(
                        filtered_signal[min_before] - filtered_signal[peak]) < 3:
                    continue
                segment = filtered_signal[min_before:min_after]
                temp = np.where(segment > np.percentile(segment, 75))[0]
                min_x = temp[0]
                max_x = temp[-1]
                tempx = []
                tempy = []
                for iii in range(min_before + min_x, min_before + max_x):
                    tempx.append(iii)
                    tempy.append(1)
                nucleusome_loc.append(int(np.mean(tempx)))
                peak_value.append(filtered_signal[peak])
                peak_loc.append(peak)
                bottom_value.append(filtered_signal[min_before])
                bottom_loc.append(min_before)
        dream_loc = [nucleusome_loc[0]]
        w_nucleu_loc.append(1)
        nucleu_loc.append(1)
        nucleusome_loc_index = 1
        while True:
            if nucleusome_loc_index >= len(nucleusome_loc) or np.mean(dream_loc[-1]) > np.mean(nucleusome_loc[-1]):
                break
            tempx = nucleusome_loc[nucleusome_loc_index]
            last_x = dream_loc[-1]
            now_x = tempx
            if now_x - last_x < 334:
                dream_loc.append(now_x)
                nucleusome_loc_index += 1
                nucleu_loc.append(1)
                w_nucleu_loc.append(1)
            else:
                now_x = last_x + 168
                dream_loc.append(now_x)
                nucleu_loc.append(0)
                w_nucleu_loc.append(0)
        while True:
            now_x = dream_loc[0]
            if now_x - 0 > 178:
                now_x = now_x - 168
                dream_loc.insert(0, now_x)
                w_nucleu_loc.insert(1, 0)
                nucleu_loc.insert(0, 0)
            else:
                break
        while True:
            now_x = dream_loc[-1]
            if len(wps_arr) - now_x > 178:
                now_x = now_x + 168
                dream_loc.append(now_x)
                w_nucleu_loc.append(0)
                nucleu_loc.append(0)
            else:
                break
        for index,_ in enumerate(dream_loc):
            if index+1 < len(dream_loc):
                w_peak_value.append(filtered_signal[dream_loc[index]])
                w_bottom_value.append(filtered_signal[int((dream_loc[index]+dream_loc[index+1])/2)])
            else:
                if int(dream_loc[index]+80) < len(dream_loc):
                    w_peak_value.append(filtered_signal[dream_loc[index]])
                    w_bottom_value.append(filtered_signal[int(dream_loc[index]+80)])
                else:
                    w_peak_value.append(filtered_signal[dream_loc[index]])
                    w_bottom_value.append(filtered_signal[-1])
            w_peak_bottom_value.append(w_peak_value[-1] - w_bottom_value[-1])
            peak_bottom_value.append(w_peak_value[-1] - w_bottom_value[-1])
        nucleusome_index= 0
        next_loc = nucleusome_loc[nucleusome_index]
        now_loc = dream_loc[0]
        for next_nucleu_exist in nucleu_loc[1:]:
            if next_nucleu_exist == 0:
                w_peak_distance.append(next_loc-now_loc)
                peak_distance.append(next_loc-now_loc)
            else:
                w_peak_distance.append(next_loc-now_loc)
                peak_distance.append(next_loc - now_loc)
                now_loc = nucleusome_loc[nucleusome_index]
                nucleusome_index+=1
                if nucleusome_index >= len(nucleusome_loc):
                    next_loc = nucleusome_loc[-1]
                else:
                    next_loc = nucleusome_loc[nucleusome_index]
        peak_bottom_value = [abs(x) for x in peak_bottom_value]
        peak_distance = [abs(x) for x in peak_distance]
        for i in range(len(nucleu_loc)):
            reference_nucleu_loc.append(1)
        for i in range(len(peak_distance)):
            reference_peak_distance.append(np.mean(peak_distance))
        for i in range(len(peak_bottom_value)):
            reference_peak_bottom_value.append(np.mean(peak_bottom_value))
        for loc1,loc2 in zip(peak_loc,bottom_loc):
            high = filtered_signal[loc1] - filtered_signal[loc2]
            k1 = high/(loc1-loc2)
            k1_list.append(abs(k1))
        for loc1, loc2 in zip(peak_loc[:-2], bottom_loc[1:]):
            high = filtered_signal[loc1] - filtered_signal[loc2]
            k2 = high / (loc1 - loc2)
            k2_list.append(abs(k2))
        for k in k1_list:
            reference_k1_list.append(np.mean(k1_list))
        for k in k2_list:
            reference_k2_list.append(np.mean(k2_list))
        if len(w_nucleu_loc[1:]) < 45:
            continue
        nucleu_loc = w_nucleu_loc[6:46]
        reference_nucleu_loc = []
        for i in range(len(nucleu_loc)):
            reference_nucleu_loc.append(1)
        if len(w_peak_bottom_value[1:]) < 45:
            continue
        peak_bottom_value = w_peak_bottom_value[6:46]
        peak_bottom_value = [abs(i) for i in peak_bottom_value]
        peak_bottom_value_mean = np.mean(np.array(peak_bottom_value))
        reference_peak_bottom_value = []
        for i in range(len(peak_bottom_value)):
            reference_peak_bottom_value.append(peak_bottom_value_mean)
        if len(w_peak_distance[1:]) < 45:
            continue
        peak_distance = w_peak_distance[6:46]
        peak_distance = [abs(i) for i in peak_distance]
        peak_distance_mean = np.mean(np.array(peak_distance))
        reference_peak_distance = []
        for i in range(len(peak_distance)):
            reference_peak_distance.append(peak_distance_mean)
        js_distance_nucleu_loc = js_distance_from_samples(nucleu_loc, reference_nucleu_loc)
        js_distance_peak_distance_value = js_distance_from_samples(peak_distance, reference_peak_distance)
        js_distance_peak_bottom_distance_value = js_distance_from_samples(peak_bottom_value, reference_peak_bottom_value)
        js_k1 = js_distance_from_samples(k1_list, reference_k1_list)
        js_k2 = js_distance_from_samples(k2_list, reference_k2_list)
        js_distance = js_distance_nucleu_loc + js_distance_peak_distance_value + js_distance_peak_bottom_distance_value + js_k1 + js_k2
        js_distance_high.append(js_distance)
        w_js_distance.append(js_distance_nucleu_loc)
        w_js_distance.append(js_distance_peak_distance_value)
        w_js_distance.append(js_distance_peak_bottom_distance_value)
        w_js_distance.append(js_k1)
        w_js_distance.append(js_k2)
        w1.writerow(w_js_distance)
        w2.writerow(w_peak_value)
        w3.writerow(w_bottom_value)
        w4.writerow(w_nucleu_loc)
        w5.writerow(w_peak_bottom_value)
        w6.writerow(w_peak_distance)
    except:
        continue
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()
