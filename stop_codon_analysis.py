#!/usr/bin/python

########################################################################
# 9 Jan 2015
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import argparse
import os, re, sys
import subprocess
from collections import defaultdict
import pybedtools 

def read_dexseq():
	dif1 = defaultdict(list)
	dif2 = defaultdict(list)
	with open("/raid/home/patrick/72_roberto_splicing/analysis/dexseq/75K.DIF_vs_CTR.DIF/75K.DIF_vs_CTR.DIF_dexseq.tsv") as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if len(word) > 5:
				if word[5] == "NA":
					pass
				else:
					pval = float(word[5])
					if pval <= 0.05:
						exon = word[0].split(":")
						dif1[word[1]].append("{},{}:{}:{}".format(exon[1], word[8], word[9], word[10]))
	with open("/raid/home/patrick/72_roberto_splicing/analysis/dexseq/75K.UND_vs_CTR.UND/75K.UND_vs_CTR.UND_dexseq.tsv") as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if len(word) > 5:
				if word[5] == "NA":
					pass
				else:
					pval = float(word[5])
					if pval <= 0.05:
						exon = word[0].split(":")
						dif2[word[1]].append("{},{}:{}-{}".format(exon[1], word[8], word[9], word[10]))
	return dif1, dif2

def get_fasta(dif1, dif2):
	fa = "/home/patrick/Reference_Genomes/hg19/UCSC/Chrom_fa/ucsc_hg19.fa"
	bed = []
	for key in sorted(dif1):
		for exon in dif1[key]:
			word = exon.split(",")
			pos = word[1].split(":")
			print "{}\t{}\t{}\n".format(pos[0], pos[1], pos[2]),

	#b = pybedtools.BedTool(bed)
	#print b
	#fa1 = b.sequence(fi=fa) 
	#print open(fa1.seqfn).read()
#	print fa1

def main():
	dif1, dif2 = read_dexseq()
	get_fasta(dif1, dif2)
	


main()