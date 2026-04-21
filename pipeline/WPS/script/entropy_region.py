#%%
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.signal import find_peaks
from scipy.stats import entropy
import argparse
parse = argparse.ArgumentParser()
parse.add_argument('-i',required=True)
parse.add_argument('-s',required=True)
args = parse.parse_args()
input_path = args.i
output_path = args.s
from tqdm import tqdm
# input_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/plasma_100/data/HC_bt31/wps.pickle'
# output_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/plasma_100/data/HC_bt31/entropy.pickle'

with open(input_path,'rb') as f:
    wps_data = pickle.load(f)
#%%
all_diffe = []
entropy_value_dict = {}
for ensg in tqdm(list(wps_data.keys())):
    wps_arr = []
    wps_pos = []
    for pos in wps_data[ensg]:
        wps_arr.append(wps_data[ensg][pos]['wps'])
        wps_pos.append(pos)
    fs = 1000  #采样率
    cutoff_frequency = 8  #截止频率
    order = 5  #滤波器阶数
    b, a = signal.butter(order, cutoff_frequency, fs=fs, btype='low')
    filtered_signal = signal.filtfilt(b, a, wps_arr)
    filtered_signal = filtered_signal - np.mean(filtered_signal)
    peaks, _ = find_peaks(filtered_signal, height=0)
    peaks = peaks.tolist()
    nucleusome_loc = []
    for i in peaks:
        left_min_index = i
        right_min_index = i
        for j in range(100):
            if i+j+1 < len(filtered_signal):
                if filtered_signal[i+j] > filtered_signal[i+j+1]:
                    right_min_index=i+j+1
                else:
                    break
            else:
                break
        for j in range(100):
            if i-j-1 > 0:
                if filtered_signal[i-j] < filtered_signal[i-j-1]:
                    left_min_index=i-j-1
                else:
                    break
            else:
                break
        segment = filtered_signal[left_min_index:right_min_index]
        seg_index = np.where(segment>0)
        if len(seg_index[0]) > 3:
            temp = np.where(segment>np.median(segment))[0]
            min_x = temp[0]
            max_x = temp[-1]
            tempx = []
            tempy = []
            for iii in range(left_min_index+min_x,left_min_index+max_x):
                tempx.append(iii)
                tempy.append(1)
            nucleusome_loc.append(np.mean(tempx))
    diff = np.diff(np.array(nucleusome_loc))
    # diff = diff//30
    # all_diffe.extend(diff)
    entropy_value = entropy(diff, base=2)  # 使用以2为底的对数计算熵
    entropy_value_dict[ensg] = entropy_value
with open(output_path, "wb") as f:
    pickle.dump(entropy_value_dict, f)
#%%
# plt.hist(diff,500)
# plt.title('30 region')
# plt.show()