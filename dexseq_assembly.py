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
import pprint
from BCBio import GFF

def read_dexseq():
	dif1 = defaultdict(list)
	dif2 = defaultdict(list)
	with open("/home/patrick/72_roberto_splicing/analysis/dexseq/75K.DIF_vs_CTR.DIF/75K.DIF_vs_CTR.DIF.tsv") as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if len(word) > 5:
				if word[6] == "NA":
					pass
				else:
					pval = float(word[6])
					if pval <= 0.05:
						exon = word[0].split(":")
						dif1[word[1]].append(exon[1])
	with open("/home/patrick/72_roberto_splicing/analysis/dexseq/75K.UND_vs_CTR.UND/75K.UND_vs_CTR.UND.tsv") as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if len(word) > 5:
				if word[6] == "NA":
					pass
				else:
					pval = float(word[6])
					if pval <= 0.05:
						exon = word[0].split(":")
						dif2[word[1]].append(exon[1])
	return dif1, dif2

def read_gff():
	in_file = "/home/patrick/72_roberto_splicing/analysis/dexseq/Homo_sapiens.GRCh37.74.ucsc_dexseq.gff"
	in_handle = open(in_file)
	for rec in GFF.parse(in_handle):
		print rec
	in_handle.close()

def main():
	read_gff()

main()