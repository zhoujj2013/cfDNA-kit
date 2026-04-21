#%%
import pysam
import os
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import csv

import argparse
import os



parser = argparse.ArgumentParser(description='处理片段分布图')
parser.add_argument('--sample_file', type=str, required=True,
                    help='样本列表文件路径 (.lst)')
parser.add_argument('--save_path', type=str, required=True,
                    help='输出图片保存路径 (例如 .png)')
parser.add_argument('--start', type=int, default=10,
                    help='起始范围 (默认: 10)')
parser.add_argument('--end', type=int, default=250,
                    help='结束范围 (默认: 250)')

args = parser.parse_args()


# 确保保存目录存在
os.makedirs(os.path.dirname(args.save_path), exist_ok=True)
sample_file = args.sample_file
save_path = args.save_path
start = args.start
end = args.end

# sample_file = '/mnt/dfc_data2/project/linyusen/project/31_cfdna_wps/pso_healthy_2026/analysis/samples.lst'
# save_path = '/mnt/dfc_data2/project/linyusen/project/31_cfdna_wps/ppt_figure/3/figure/fragment.png'
# start = 10
# end = 250


def require_fragment_length(bam_list,range):
    length_dict = {}
    for bam_path in bam_list:
        file = os.path.basename(bam_path)
        print(file)
        bam_name = file.split('.')[0]
        length_dict[bam_name] = []
        sf = pysam.AlignmentFile(bam_path)
        for index,reads in enumerate(sf):
            if index == 10000000:
                break
            if abs(reads.isize) > range[1]:
                continue
            if abs(reads.isize) < range[0]:
                continue
            fragment_length = abs(reads.isize)
            length_dict[bam_name].append(fragment_length)
    return length_dict


bam_list = {}
f = open(sample_file)
r = csv.reader(f)
for i in r:
    group,bam_file = i
    if group not in bam_list:
        bam_list[group] = []
    bam_list[group].append(bam_file)
range = [start,end]


fragment_length = {}
for group_name in bam_list:
    fragment_length[group_name] = require_fragment_length(bam_list[group_name],range)

colors = ['red','blue']
plt.figure(figsize=[8,8],dpi=300)
# 示例数据
for index,group_name in enumerate(fragment_length):
    for bam_name in fragment_length[group_name]:
        data = fragment_length[group_name][bam_name]
        # 绘制概率密度分布图
        sns.kdeplot(data, linewidth=1,label=group_name,color=colors[index])

bwith = 2
# 添加标题和标签
plt.title('Probability Density Distribution')
plt.xlabel('Value')
plt.ylabel('Density')
ax = plt.gca()#获取边框
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
ax.spines['bottom'].set_linewidth(bwith)
plt.legend()
plt.tick_params(axis='both', which='major', width=bwith, length=bwith*3,labelsize=12)
plt.savefig(save_path)
plt.show()
