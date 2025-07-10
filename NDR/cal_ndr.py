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

# 


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
		flaglen = abs(int(x.template_length))
		if x.flag == 99 and x.next_reference_id == x.reference_id:
			# we will remove the flag size > 2k, which will not be seq in MGI
			if flaglen < 2000:
				xs = x.reference_start
				xe = x.reference_start + abs(int(x.template_length))
				for p in range(xs,xe,1):
					rpos = p - tss
					if rpos > -1000 and rpos < 1000:
						pos.append(rpos)
				i = i + 1
				#depth=depth+abs(int(x.template_length))
		elif x.flag == 83 and x.next_reference_id == x.reference_id:
			if flaglen < 2000:
				xs = x.reference_start
				xe = x.reference_start - abs(int(x.template_length))
				for p in range(xe,xs,1):
					rpos = p - tss
					if rpos > -1000 and rpos < 1000:
						pos.append(rpos)
				i = i + 1
				#depth=depth+abs(int(x.template_length))
	
	# calculate flagment count in a specific region
	pos.append(0)

	vect_pos = np.histogram(pos, binedges)

	## get upstream mean depth	
	utss=int(c[1])-2000
	us=int(c[1])-3000
	ue=int(c[2])-1000
	
	upos = []
	i = 0
	#depth=0
	iter = sf.fetch(chrid, us, ue)
	for x in iter:
		# 99 read 1 forward, properly paired
		# 83 read 1 reverse, properly paired
		flaglen = abs(int(x.template_length))
		if x.flag == 99 and x.next_reference_id == x.reference_id:
			# we will remove the flag size > 2k, which will not be seq in MGI
			if flaglen < 2000:
				xs = x.reference_start
				xe = x.reference_start + abs(int(x.template_length))
				for p in range(xs,xe,1):
					rpos = p - utss
					if rpos > -1000 and rpos < 1000:
						upos.append(rpos)
				i = i + 1
				#depth=depth+abs(int(x.template_length))
		elif x.flag == 83 and x.next_reference_id == x.reference_id:
			if flaglen < 2000:
				xs = x.reference_start
				xe = x.reference_start - abs(int(x.template_length))
				for p in range(xe,xs,1):
					rpos = p - utss
					if rpos > -1000 and rpos < 1000:
						upos.append(rpos)
				i = i + 1
				#depth=depth+abs(int(x.template_length))
	
	upos.append(0)
	vect_upos = np.histogram(upos, binedges)
	upMean = np.mean(vect_upos[0])

	## get downstream signal	
	dtss=int(c[1])+2000
	ds=int(c[1])+1000
	de=int(c[2])+3000
	
	dpos = []
	i = 0
	#depth=0
	iter = sf.fetch(chrid, ds, de)
	for x in iter:
		# 99 read 1 forward, properly paired
		# 83 read 1 reverse, properly paired
		flaglen = abs(int(x.template_length))
		if x.flag == 99 and x.next_reference_id == x.reference_id:
			# we will remove the flag size > 2k, which will not be seq in MGI
			#if x.template_length < 2000:
			if flaglen < 2000:
				xs = x.reference_start
				xe = x.reference_start + abs(int(x.template_length))
				for p in range(xs,xe,1):
					rpos = p - dtss
					if rpos > -1000 and rpos < 1000:
						dpos.append(rpos)
				i = i + 1
				#depth=depth+abs(int(x.template_length))
		elif x.flag == 83 and x.next_reference_id == x.reference_id:
			if flaglen < 2000:
				xs = x.reference_start
				xe = x.reference_start - abs(int(x.template_length))
				for p in range(xe,xs,1):
					rpos = p - dtss
					if rpos > -1000 and rpos < 1000:
						dpos.append(rpos)
				i = i + 1
				#depth=depth+abs(int(x.template_length))
	
	dpos.append(0)
	vect_dpos = np.histogram(dpos, binedges)

	# final mean value
	downMean = np.mean(vect_dpos[0])
	
	updownMean = (upMean + downMean)/2;
	
	vectt_pos = [float(depth/updownMean) for depth in vect_pos[0]]
	
	# print the result
	print("\t".join(c), end="")
	print("\t", end="")
	
	# print positve strand signal
	#print(",".join([str(i) for i in vectt_ppos]), end="")
	#print("\t", end="")

	# print negative strand signal	
	print(",".join([str(i) for i in vectt_pos]))
bedfh.close()
