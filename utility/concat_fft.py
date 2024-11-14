import os
import gzip
import tqdm
import argparse
import pickle

parse = argparse.ArgumentParser()
parse.add_argument('-i',required=True,type=str)
parse.add_argument('-s',required=True,type=str)
args = parse.parse_args()
input_path = args.i
output_path = args.s

# input_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/ranking_test/data/HC_zhou2_04/fft'
# output_path = '/mnt/dfc_data2/project/linyusen/database/46_cfdna/newdata/ranking_test/data/HC_zhou2_04/fft.pickle'
data = {}
for file in tqdm.tqdm(os.listdir(input_path)):

    gene = file.split('.')[0].split('_')[1]
    data[gene] = {}
    with gzip.open(os.path.join(input_path, file), 'rb') as f:
        content = f.read()
    content_str = content.decode('utf-8').split('\n')
    for line in content_str[1:-1]:
        line = line.split('\t')
        fre = int(line[0])
        intensity = float(line[-1])
        data[gene][fre] = intensity
#%%

with open(output_path,'wb') as pickle_file:
    pickle.dump(data, pickle_file)