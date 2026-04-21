#!/usr/bin/env python3

import sys, os
from optparse import OptionParser
import gzip
import pysam
import random
from collections import defaultdict
from bx.intervals.intersection import Intersecter, Interval

def isSoftClipped(cigar):
    for (op, count) in cigar:
        if op in (4, 5, 6):
            return True
    return False

def aln_length(cigarlist):
    tlength = 0
    for operation, length in cigarlist:
        if operation in (0, 2, 3) or operation >= 6:
            tlength += length
    return tlength

parser = OptionParser()
parser.add_option("-i","--input", dest="input", default="transcriptAnno.tsv")
parser.add_option("-l","--lengthSR", dest="lengthSR", default=76, type="int")
parser.add_option("-m","--merged", dest="merged", default=False, action="store_true")
parser.add_option("-t","--trimmed", dest="trimmed", default=False, action="store_true")
parser.add_option("-w","--protection", dest="protection", default=120, type="int")
parser.add_option("-o","--outfile", dest="outfile", default='block_%s.tsv.gz')
parser.add_option("-e","--empty", dest="empty", default=False, action="store_true")
parser.add_option("--minInsert", dest="minInsSize", default=-1, type="int")
parser.add_option("--maxInsert", dest="maxInsSize", default=-1, type="int")
parser.add_option("--max_length", dest="max_length", default=1000, type="int")
parser.add_option("--downsample", dest="downsample", default=None, type="float")
parser.add_option("-v","--verbose", dest="verbose", default=False, action="store_true")
(options, args) = parser.parse_args()

minInsSize = maxInsSize = None
if options.minInsSize > 0 and options.maxInsSize > 0 and options.minInsSize < options.maxInsSize:
    minInsSize = options.minInsSize
    maxInsSize = options.maxInsSize

protection = options.protection // 2
validChroms = []
for i in range(1,23):
    validChroms.append(str(i))
validChroms.extend(["X","Y"])


if os.path.exists(options.input):
    with open(options.input) as infile:
        for line in infile:
            cid, chrom, start, end, strand = line.split()
            if chrom not in validChroms:
                continue

            regionStart, regionEnd = int(start), int(end)
            if regionStart < 1:
                continue

            posRange = defaultdict(lambda:[0,0])
            filteredReads = Intersecter()

            for bamfile in args:
                bamfile = bamfile.strip("'")
                if not os.path.exists(bamfile):
                    continue

                bam = pysam.AlignmentFile(bamfile, "rb")
                prefix = "chr" if any(r.startswith("chr") for r in bam.references) else ""

                for read in bam.fetch(prefix+chrom,
                                      regionStart-protection-1,
                                      regionEnd+protection+1):

                    if read.is_duplicate or read.is_qcfail or read.is_unmapped:
                        continue
                    if read.cigar is None or isSoftClipped(read.cigar):
                        continue

                    if read.is_paired:
                        if read.mate_is_unmapped:
                            continue
                        if read.next_reference_id != read.reference_id:
                            continue
                        if read.is_read1 or (
                            read.is_read2 and
                            read.next_reference_start + read.query_length < regionStart-protection-1
                        ):
                            if read.template_length == 0:
                                continue
                            if options.downsample and random.random() >= options.downsample:
                                continue

                            rstart = min(read.reference_start,
                                         read.next_reference_start) + 1
                            lseq = abs(read.template_length)
                            rend = rstart + lseq - 1

                            if minInsSize and not (minInsSize <= lseq <= maxInsSize):
                                continue

                            filteredReads.add_interval(Interval(rstart, rend))
                            for i in range(rstart, rend+1):
                                if regionStart <= i <= regionEnd:
                                    posRange[i][0] += 1
                            if regionStart <= rstart <= regionEnd:
                                posRange[rstart][1] += 1
                            if regionStart <= rend <= regionEnd:
                                posRange[rend][1] += 1

                    else:
                        if options.downsample and random.random() >= options.downsample:
                            continue

                        rstart = read.reference_start + 1
                        lseq = aln_length(read.cigar)
                        rend = rstart + lseq - 1

                        if minInsSize and not (minInsSize <= lseq <= maxInsSize):
                            continue

                        filteredReads.add_interval(Interval(rstart, rend))
                        for i in range(rstart, rend+1):
                            if regionStart <= i <= regionEnd:
                                posRange[i][0] += 1

                        if ((options.merged or read.query_name.startswith('M_')) or
                            ((options.trimmed or read.query_name.startswith('T_')) and
                             read.query_length <= options.lengthSR-10)):
                            if regionStart <= rstart <= regionEnd:
                                posRange[rstart][1] += 1
                            if regionStart <= rend <= regionEnd:
                                posRange[rend][1] += 1
                        elif read.is_reverse:
                            if regionStart <= rend <= regionEnd:
                                posRange[rend][1] += 1
                        else:
                            if regionStart <= rstart <= regionEnd:
                                posRange[rstart][1] += 1

            filename = options.outfile % cid
            cov_sites = 0

            with gzip.open(filename, 'wt') as outfile:
                outLines = []
                for pos in range(regionStart, regionEnd+1):
                    rstart, rend = pos-protection, pos+protection
                    gcount = bcount = 0
                    for read in filteredReads.find(rstart, rend):
                        if read.start > rstart or read.end < rend:
                            bcount += 1
                        else:
                            gcount += 1

                    covCount, startCount = posRange[pos]
                    cov_sites += covCount
                    outLines.append(
                        f"{chrom}\t{pos}\t{covCount}\t{startCount}\t{gcount-bcount}\n"
                    )

                if strand == "-":
                    outLines.reverse()

                outfile.writelines(outLines)

            if cov_sites == 0 and not options.empty:
                os.remove(filename)
