#%%
import os
from scipy.signal import find_peaks
from scipy import signal
import numpy as np
import pysam
import argparse
from tqdm import tqdm
import csv
from scipy.spatial.distance import jensenshannon


parser = argparse.ArgumentParser(description="Process samples and BED file to generate RPM results")
parser.add_argument(
    "--sample_file",
    type=str,
    required=True,
    help="Path to the sample list file (one sample per line)"
)
parser.add_argument(
    "--bed_file",
    type=str,
    required=True,
    help="Path to the BED file with genomic regions"
)
parser.add_argument(
    "--save_dir",
    type=str,
    required=True,
    help="Directory to save the output results"
)
args = parser.parse_args()

sample_file = args.sample_file
bed_file = args.bed_file
save_dir = args.save_dir


# sample_file = '/mnt/dfc_data2/project/linyusen/project/31_cfdna_wps/project/concate_region/samples.lst'
# bed_file = '/mnt/dfc_data2/project/zhoujj/project/35cfdna/04PanCancer/ref/gencode.v45.annotation.tss.b5k.bed'
# save_dir = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/concate_region/pfe'
os.makedirs(save_dir, exist_ok=True)
#%%
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
def cal_nucleu_loc(wps_arr):
    nucleu_loc = []
    wps_arr = np.array(wps_arr)
    b, a = signal.butter(5, 8, fs=1000, btype='low')
    filtered_signal = signal.filtfilt(b, a, wps_arr)
    filtered_signal = filtered_signal - np.mean(filtered_signal)
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

    dream_loc = [nucleusome_loc[0]]
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
        else:
            now_x = last_x + 168
            dream_loc.append(now_x)
            nucleu_loc.append(0)
    while True:
        now_x = dream_loc[0]
        if now_x - 0 > 178:
            now_x = now_x - 168
            dream_loc.insert(0, now_x)
            nucleu_loc.insert(0, 0)
        else:
            break
    while True:
        now_x = dream_loc[-1]
        if len(wps_arr) - now_x > 178:
            now_x = now_x + 168
            dream_loc.append(now_x)
            nucleu_loc.append(0)
        else:
            break
    return nucleu_loc

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
region_data = region_data[:100]

f = open(sample_file)
r = csv.reader(f,delimiter=' ')
sample_dict = {}
for line in r:
    group = line[0]
    name = line[1]
    path = line[2]
    if group not in sample_dict:
        sample_dict[group] = {}
    sample_dict[group][name] = path
f.close()

window = 1
js_dict = {}
for group in sample_dict:
    js_dict[group] = {}
    for sample_name in sample_dict[group]:

        bam_file = sample_dict[group][sample_name]
        nucleu_list = []
        reference_nucleu_list = []
        for chr, start, end, gene in tqdm(region_data,desc=sample_name):
            wps_arr, wps_arr_x = wps_signal(bam_file, 120, chr, start, end)
            nucleu_loc = cal_nucleu_loc(wps_arr)
            length = len(nucleu_loc)
            median_loc = int(length / 2)
            median_nucleu = nucleu_loc[median_loc - int(window / 2):median_loc + int(window / 2 + 0.5)]
            nucleu_list.extend(median_nucleu)
        for i in nucleu_list:
            reference_nucleu_list.append(1)
        js_distance = jensenshannon(nucleu_list, reference_nucleu_list)

        js_dict[group][sample_name] = js_distance
#%%
for group in js_dict:
    f = open(os.path.join(save_dir, group + '.js.tsv'),'w')
    w = csv.writer(f,delimiter='\t')
    for sample_name in js_dict[group]:
        w.writerow([sample_name, js_dict[group][sample_name]])
    f.close()


