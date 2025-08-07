import pysam
import numpy as np
from scipy.stats import entropy

def extract_fragment_lengths(bamfile, chrom, start, end):
    """Extract fragment lengths from properly paired reads within region."""
    lengths = []
    for read in bamfile.fetch(chrom, start, end):
        if read.flag in (99, 83) and read.reference_id == read.next_reference_id:
            if abs(read.template_length) < 2000:
                lengths.append(abs(read.template_length))
    return lengths
def split_region(start, end, window_size):
    """
    将给定区域 [start, end) 按 window_size 划分多个子区域。

    返回：[(start1, end1), (start2, end2), ...]
    """
    regions = []
    for i in range(start, end, window_size):
        sub_start = i
        sub_end = min(i + window_size, end)
        regions.append([sub_start, sub_end])
    return regions
def compute_entropy(lengths, bins):
    """Compute histogram and entropy of fragment length distribution."""
    filtered_lengths = [l for l in lengths if 50 <= l <= 400]
    if len(filtered_lengths) < 20:
        return None  # too few fragments
    hist, _ = np.histogram(filtered_lengths, bins=bins)
    probs = hist / np.sum(hist)
    return entropy(probs)


def pfe_signal(sf:pysam.AlignmentFile, windows, chrid, start, end):
    # bamfile = pysam.AlignmentFile(bam_file, "rb")
    regions = split_region(start,end,windows)
    bins = list(range(50, 400, 50))
    entropy_vals = []
    for sub_start, sub_end in regions:
        lengths = extract_fragment_lengths(sf, chrid, sub_start, sub_end)
        entropy_val = compute_entropy(lengths,bins)
        entropy_vals.append(entropy_val)
    return entropy_vals
