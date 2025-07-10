import os, sys
import re

import pysam
from scipy.stats import entropy
import numpy as np

def usage():
	print('\nCalculate pfe value (entropy) for ROI in cfDNA ngs\n')
	print('Author: zhoujj2013@gmail.com 6/2023')
	print('Usage: python '+sys.argv[0]+' xxx.bam xxx.bed\n')
	print('Note:')
	print('xxxx.bed contain regions (for instance, arround TSS, SE, TE etc.) in bed format.')
	print('Require: > python3')
	print('')
	sys.exit(2)

# check args
if len(sys.argv) < 2:
	usage()

bamf = sys.argv[1]
bedf = sys.argv[2]

disStart = int(sys.argv[3])
disEnd = int(sys.argv[4])

sf = pysam.AlignmentFile(bamf, "rb")
binedges = range(50, 400, 50) # bin setting is very important parameters----jiajian

bedfh = open(bedf,'r')
while True:
	l = bedfh.readline()
	
	if len(l) == 0:
        	break
	
	l = l.strip('\n')

	#print(l)
	c = l.split('\t')
	
	strand = c[5]
	chrid=c[0]
	s=int(c[1])
	e=int(c[2])
	ss = 0
	ee = 0

	if strand == "+":
		ss = s + disStart
		ee = e + disEnd
	elif strand == "-":
		ss = s - disEnd
		ee = e - disStart

	# maybe we do not need this
	#ss = ss - 500
	#ee = ee - 500
	#print(chrid + "\t" + str(ss) + "\t" + str(ee) + "\t" + c[3] + "\t" + c[5])

	#print(chrid + "##" + str(s) + "##" + str(e))
	fsize = {}
	iter = sf.fetch(chrid, ss, ee)
	#print(len(iter))
	for x in iter:
		#print(str(x.next_reference_id) + "\t" + str(x.reference_id))
		# only keep properly paired reads, edited by jiajian, 2023.06
		#if (x.flag == 99 or x.flag == 97 or x.flag == 147 or x.flag == 145 or x.flag == 163 or x.flag == 161 or x.flag == 83 or x.flag == 81) and x.next_reference_id == x.reference_id:
		if (x.flag == 99 or x.flag == 83) and x.next_reference_id == x.reference_id:
			# we will remove the flag size > 2k, which will not be seq in MGI
			if x.template_length < 2000:
				fsize[x.query_name] = abs(x.template_length)
			#print(x.query_name + "\t" + str(abs(x.template_length)))
			#print(str(x.template_length))
	# calculate density of a specific region, you can output histogram to check the distribution.
	vvlaueArr = [(i <= 400 and i >= 50) for i in fsize.values()]
	if len(vvlaueArr) <= 20: # min fragments in the windows
		continue
	#print("#".join([str(i) for i in fsize.values()]), end="\n")
	lenDensity = np.histogram([float(i) for i in fsize.values()], binedges, density=False)
	count_sum = np.sum(lenDensity[0])
	lenDensityCount = [i/count_sum for i in lenDensity[0]]
	#print("#".join([str(i) for i in lenDensityCount]), end="\n")
	
	# calculate the entropy of a distribution for given probability values.
	tssEntropy = entropy(lenDensityCount)
	# 2023/6/6, zhoujj and sunhang edit
	# it require high sequencing deep for a specific region.

	# print the result
	print("\t".join(c), end="")
	print("\t", end="")
	#print(" ".join([str(i) for i in lenDensity[0]]), end="")
	#print("\t", end="")
	print(tssEntropy)
bedfh.close()
