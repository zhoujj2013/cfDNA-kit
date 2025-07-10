import os, sys
import re

import pysam
from scipy.stats import entropy
import numpy as np

def usage():
	print('\nCalculate the nucleosome-depleted region signal for TSS (+/-1k) using cfDNA NGS\n')
	print('Author: zhoujj2013@gmail.com 5/12/2023\n')
	print('Usage: python '+sys.argv[0]+' xxx.bam xxx.tss.bed  > ndr.result 2>total.flagment.count')
	print('')
	sys.exit(2)

# check args
if len(sys.argv) < 2:
	usage()

bamf = sys.argv[1]
bedf = sys.argv[2]



# get properly paired mapped flagment
'''
result = pysam.flagstat("-@","16", bamf)
pattern = re.compile(r'(\d+) \+ (\d+) properly paired', re.M|re.I)
m = pattern.search(result)
properly_paired_mapped = m.group(1)

print("Total fragments: " + properly_paired_mapped, file=sys.stderr)
'''

#print(properly_paired_mapped)

sf = pysam.AlignmentFile(bamf, "rb")

binedges = range(-1000, 1000, 5)

bedfh = open(bedf,'r')
while True:
	l = bedfh.readline()
	
	if len(l) == 0:
        	break
	
	l = l.strip('\n')

	#print(l)
	c = l.split('\t')
	chrid=c[0]
	tss=int(c[1])
	s=int(c[1])-1500
	e=int(c[2])+1500
	
	pos = []
	i = 0
	#depth=0
	iter = sf.fetch(chrid, s, e)
	for x in iter:
		# 99 read 1 forward, properly paired
		# 83 read 1 reverse, properly paired
		# for single end sequencing
		if x.flag == 0:
			# we will remove the flag size > 2k, which will not be seq in MGI
			xs = x.reference_start
			xe = x.reference_start + 35
			for p in range(xs,xe,1):
				rpos = p - tss
				if rpos > -1000 and rpos < 1000:
					pos.append(rpos)
			i = i + 1
			#depth=depth+abs(int(x.template_length))
		elif x.flag == 16:
			xs = x.reference_start
			xe = x.reference_start - 35
			for p in range(xe,xs,1):
				rpos = p - tss
				if rpos > -1000 and rpos < 1000:
					pos.append(rpos)
			i = i + 1
			#depth=depth+abs(int(x.template_length))
	
	# calculate flagment count in a specific region
	pos.append(0)

	vect_pos = np.histogram(pos, binedges)

	vmin = min(vect_pos[0])
	vmax = max(vect_pos[0])

	vectt_pos = [float((depth-vmin)/(vmax-vmin)) for depth in vect_pos[0]]
	
	# print the result
	print("\t".join(c), end="")
	print("\t", end="")
	
	# print positve strand signal
	#print(",".join([str(i) for i in vectt_ppos]), end="")
	#print("\t", end="")

	# print negative strand signal	
	print(",".join([str(i) for i in vectt_pos]))
bedfh.close()
