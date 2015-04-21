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
import ConfigParser
from multiprocessing import Pool
import itertools

def index(fa, index):
	#First download genes and gene predictions ensembl table from UCSC and save as hg19_transcriptome.fa 
	command = "salmon index -t {} -i {}".format(fa, index)
	subprocess.call(command.split())

def quant(fq1, threads, outdir, lib, fq2=None):
	fa = "/home/patrick/Reference_Genomes/hg19/UCSC/hg19_transcriptome.fa"
	index = "/home/patrick/Reference_Genomes/hg19/UCSC/salmon_index"
	gtf = "/home/patrick/72_roberto_splicing/annotation/Homo_sapiens.GRCh37.74.ucsc.gtf"
	#For current, use ISR
	if fq2:
		command2 = "salmon quant -p {} -i {} -l {} -1 {} -2 {} -o {}".format(threads, index, lib, fq1, fq2, outdir)
	else:
		command2 = "salmon quant -p {} -i {} -l SR -r {} -o {}".format(threads, index, fq1, outdir)
	subprocess.call(command2.split())

def combine_quant(conditions, outfile):
	data = defaultdict(list)
	output = open(outfile, "w")
	output.write("ID"),
	for bam in sorted(conditions):
		name = os.path.basename(bam)
		name = name.strip(".bam$")
		output.write("\t{}".format(name)),
		with open("{}_salmon/quant.sf".format(name)) as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				if word[0].startswith("#"):
					pass
				else:
					data[word[0]].append(word[2])
	output.write("\n"),
	for key2 in sorted(data):
		data2 = data[key2]
		output.write(key2+"\t" + "\t".join(data2) + "\n"),
	output.close()

def ConfigSectionMap(section, Config):
	dict1 = {}
	options = Config.options(section)
	for option in options:
		try:
			dict1[option] = Config.get(section, option)
			if dict1[option] == -1:
				DebugPrint("skip: %s" % option)
		except:
			print("exception on %s!" % option)
			dict1[option] = None
	return dict1

def quant_fun(args):
	return quant(*args)

def main():
	parser = argparse.ArgumentParser(description='Transcript quantification using salmon\n')
	parser.add_argument('-f', '--fastq', help='Single end fastq', required=False)
	parser.add_argument('-p', '--pair', help='Paired end fastqs. Please put them in order!', required=False, nargs='+')
	parser.add_argument('-t','--threads', help='threads, default=1', default=1, required=False)
	parser.add_argument('-l', '--lib', help='libtype for paired samples, see docs', required=False)
	parser.add_argument('-o', '--out', help='Output', required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())


	if args["fastq"]:
		quant(args["fastq"], args["threads"], args["out"], lib=None, fq2=None)
	else:
		quant(args["pair"][0], args["threads"], args["out"], args["lib"], args["pair"][1])

	#combine_quant(conditions, args["out"])

main()