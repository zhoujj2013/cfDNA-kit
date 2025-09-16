import os, sys
import re

def usage():
    print('\nCombine cufflink quantification result to matrix\n')
    print('Author: zhoujj2013@gmail.com 8/29/2016\n')
    print('Usage: python '+sys.argv[0]+' 1.gene_tracking:aa 2.gene_tracking:bb ...')
    print('')
    sys.exit(2)

# check args
if len(sys.argv) < 3:
	usage()


import pandas as pd
import scipy as sci
import numpy as np
import csv

df = pd.read_csv(sys.argv[1], sep="\t", header=None, index_col=0)

out = []
for index, row in df.iterrows():
    outrow = []
    rowl = row.to_list()
    ctl1 = rowl[0:6]
    ctl2 = rowl[6:12]
    ctl3 = rowl[3:9]
    pdac = rowl[12:]

    r1 = sci.stats.ttest_ind(ctl1, pdac)
    r2 = sci.stats.ttest_ind(ctl2, pdac)
    r3 = sci.stats.ttest_ind(ctl3, pdac)
    #print(str(r[0]) + "\t" + str(r[1]))

    outrow.append(index)
    outrow.extend([str(i) for i in r1])
    outrow.extend([str(i) for i in r2])
    outrow.extend([str(i) for i in r3])
    out.append(outrow)
    #print(index + "\t" + "\t".join([str(i) for i in zs]))


# In[94]:


with open(sys.argv[2], 'w', newline='', encoding='UTF8') as csvfile:
    writer = csv.writer(csvfile,delimiter ='\t')
    writer.writerows(out)
