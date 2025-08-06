#%%
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from scipy.fft import fft, fftfreq
from scipy import signal
import numpy as np
import pysam
import argparse
from sklearn.decomposition import NMF
import numpy as np
from sklearn.linear_model import Lasso


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
def process_region(chrid, start, end, gene,gene_list):
    if gene not in gene_list:
        return None  # 跳过不感兴趣的基因
    if start > end:
        start, end = end, start
    # 每个线程单独打开 BAM 文件（线程安全）
    wps_arr, wps_arr_x = wps_signal(bam_file, 120, chrid, start, end)
    return gene, wps_arr
def compute_band_mean_powers(wps_arr, fs=1000):
    """
    对 wps_arr 做傅里叶变换，提取 0-250Hz，每 5Hz 一组的均值功率
    :param wps_arr: 一维信号数组
    :param fs: 采样频率，默认 1000Hz（你要根据实际设置）
    :return: 每个 5Hz 频段的功率均值列表
    """
    N = len(wps_arr)
    fft_vals = fft(wps_arr)
    freqs = fftfreq(N, d=1/fs)

    # 取正频率部分
    pos_mask = freqs >= 0
    fft_vals = fft_vals[pos_mask]
    freqs = freqs[pos_mask]

    # 计算功率谱
    power = np.abs(fft_vals) ** 2

    # 提取 0-250Hz 的部分
    freq_mask = freqs <= 250
    freqs = freqs[freq_mask]
    power = power[freq_mask]

    # 分成每 5Hz 一段，计算每段均值
    band_means = []
    for start in range(0, 250, 10):
        band_mask = (freqs >= start) & (freqs < start + 10)
        if np.any(band_mask):
            mean_power = np.mean(power[band_mask])
        else:
            mean_power = 0
        band_means.append(mean_power)

    return band_means

def cal_nmf_weight(input_matrix,tissue_expression_matrix):
    n_components = 10
    # Create the NMF model with specified number of components and random initialization
    nmf = NMF(n_components=n_components, init='random', random_state=0)
    # Fit the NMF model to the input matrix and obtain the factorized matrices
    L = nmf.fit_transform(input_matrix)
    W = nmf.components_
    C = np.zeros((24, n_components))
    alpha = 0.0001
    for j in range(L.shape[1]):
        y = L[:, j]
        X = tissue_expression_matrix  # Feature matrix
        lasso = Lasso(alpha=alpha, fit_intercept=False, max_iter=200000)
        lasso.fit(X, y)
        C[:, j] = lasso.coef_
    Weight = np.dot(C, W)
    # Reconstruct the original energy matrix using the tissue expression matrix and Weight
    wps_energies_matrix_hat = tissue_expression_matrix @ Weight
    error = np.linalg.norm(input_matrix - wps_energies_matrix_hat) / np.linalg.norm(input_matrix)
    return Weight


def parse_args():
    parser = argparse.ArgumentParser(description="""
    cfDNA Tissue/Cell Contribution Estimation Module:
    This program integrates cfDNA sequencing data and single-cell gene expression matrices
    to estimate the contribution of different tissues/cell types to the cfDNA detected in blood.
    """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to the average gene expression matrix from single-cell data (e.g., Tabula Sapiens RNA)')
    parser.add_argument('--bed_file', type=str, required=True,
                        help='Path to the BED file defining target gene regions (e.g., upstream promoters)')
    parser.add_argument('--bam_file', type=str, required=True,
                        help='Path to the cfDNA aligned BAM file')
    parser.add_argument('--csv_file', type=str, required=True,
                        help='Output path for the resulting energy weight matrix (CSV format)')
    parser.add_argument('--ref_chrid', type=str, default='chr7',
                        help='Reference chromosome ID (e.g., chr7)')
    parser.add_argument('--ref_start', type=int, default=73006144,
                        help='Start coordinate of the reference region (e.g., 73006144)')
    parser.add_argument('--ref_end', type=int, default=73016144,
                        help='End coordinate of the reference region (e.g., 73016144)')
    parser.add_argument('--windows', type=int, default=120,
                        help='Sliding window size in base pairs (default: 120bp)')

    return parser.parse_args()

args = parse_args()
data_path = args.data_path
bed_file = args.bed_file
bam_file = args.bam_file
csv_file = args.csv_file
ref_chrid = args.ref_chrid
ref_start = args.ref_start
ref_end = args.ref_end
windows = args.windows

# data_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/nc_data_single_cell/Tabula_Sapiens_All_cells_averages_assays_RNA_data.txt'
# bed_file = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/cfDNA-master/expression/transcriptAnno-hg38.75.upstream.bed'
# bam_file = '/mnt/dfc_data2/project/zhoujj/project/35cfdna/00baseline/cohort2_out/PL230613017/01alignment/PL230613017.sorted.rmdup.bam'
# csv_file = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/nmf_data/PL230613017.weight.csv'
# ref_chrid = 'chr7'
# ref_start = 73006144
# ref_end = 73016144
# windows = 120




fs = 1000  # 采样率
cutoff_frequency = 100  # 截止频率
order = 10  # 滤波器阶数
b, a = signal.butter(order, cutoff_frequency, fs=fs, btype='low')
ref_wps_arr,ref_wps_arr_x = wps_signal(bam_file,windows,ref_chrid,ref_start,ref_end)
ref_wps_arr = signal.filtfilt(b, a, ref_wps_arr)
ref_wps_arr = ref_wps_arr - np.mean(ref_wps_arr)
ref_len = ref_end - ref_start
#%%
df = pd.read_csv(data_path, sep='\t')
dict_data = df.to_dict('dict')
expression_martrix = {}
for key in dict_data:
    i = key.split('#')[-1]
    expression_martrix[i] = {}
    for j in dict_data[key]:
        expression_martrix[i][j] = 0
for key in dict_data:
    i = key.split('#')[-1]
    for j in dict_data[key]:
        expression_martrix[i][j]+=dict_data[key][j]
tissues_list = sorted(expression_martrix.keys())
gene_list = sorted({f for d in expression_martrix.values() for f in d.keys()})

#%%
f = open(bed_file)
lines = f.readlines()
f.close()
tsv_data = []
for line in lines:
    tsv_data.append(line[:-1].split('\t'))

gene_wps_dict = {}
# 多线程处理
with ThreadPoolExecutor(max_workers=8) as executor:
    futures = []
    for chrid, start, end, gene, _ in tsv_data[:]:
        futures.append(executor.submit(process_region, chrid, int(start), int(end), gene,gene_list))

    for future in tqdm(as_completed(futures), total=len(futures)):
        result = future.result()
        if result:
            gene, wps_arr = result
            if len(wps_arr) != ref_len:
                continue
            gene_wps_dict[gene] = wps_arr
#%%
genes = sorted(gene_wps_dict.keys())
uniform_genes = [g for g in genes if len(gene_wps_dict[g]) == ref_len]
tissue_expression_matrix = np.zeros((len(tissues_list), len(uniform_genes)))
for i, tissue in enumerate(tissues_list):
    for j, feature in enumerate(uniform_genes):
        tissue_expression_matrix[i, j] = expression_martrix[tissue].get(feature, 0)
tissue_expression_matrix = tissue_expression_matrix.T
#%%
wps_energies_matrix = []
for gene in uniform_genes:
    wps_arr = gene_wps_dict[gene]
    wps_arr = signal.filtfilt(b, a, wps_arr)
    band_means = compute_band_mean_powers(wps_arr)  # Compute mean power values in each predefined frequency band from the filtered error signal
    wps_energies_matrix.append(band_means)
wps_energies_matrix = np.array(wps_energies_matrix)
#%%

Weight = cal_nmf_weight(wps_energies_matrix,tissue_expression_matrix)
weight_dict = {}
for i,value in enumerate(Weight):
    weight_dict[tissues_list[i]] = {}
    for j,v in enumerate(Weight[i]):
        weight_dict[tissues_list[i]][j] = v
weight_dict = pd.DataFrame(weight_dict).T
weight_dict.to_csv(csv_file)
