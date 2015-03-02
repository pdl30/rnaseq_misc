#!/usr/bin/python

########################################################################
# 19 May 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, re, os
import ConfigParser
import itertools
import argparse
from collections import defaultdict
from pyrnatools.tools import  deseq2
import pysam

def create_design_for_R(idict):
	output = open("tmp_design.txt", "w")
	output.write("sampleName\tfileName\tcondition\n"),
	for key in sorted(idict.keys()):
		name = re.sub(".bam$", "", key)
		count = re.sub(".bam$", ".count", key)
		output.write("{}\t{}\t{}\n".format(name, count, idict[key]))
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

def run_rcode(rscript, name):
	rcode = open(name, "w")
	rcode.write(rscript)
	rcode.close()
	try:
		subprocess.call(['Rscript', name])
	except:
		error("Error in running {}\n".format(name))
		error("Error: %s\n" % str(sys.exc_info()[1]))
		error( "[Exception type: %s, raised in %s:%d]\n" % ( sys.exc_info()[1].__class__.__name__, 
		os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
		traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
		sys.exit(1)

def print_norm_counts(counts_file):
	rscript =  "library(DESeq2)\n"
	rscript += "pdata <- read.table('tmp_design.txt', header=T)\n"
	#Make sure they match!
	rscript += "counts <- read.table('{}', sep='\\t', header=T, row.names=1)\n".format(counts_file)
	rscript += "rnaseq_dds <- DESeqDataSetFromMatrix(countData = counts, colData = data.frame(pdata), design = ~ condition)\n"
	rscript += "rnaseq_dds$condition <- factor(rnaseq_dds$condition, levels=unique(pdata[,3]))\n"
	rscript += "rnaseq_dds <- DESeq(rnaseq_dds)\n"
	rscript += "factors <- sizeFactors(rnaseq_dds)\n"
	rscript += "names(factors) <- rnaseq_dds$condition\n"
	rscript += "co <- counts(rnaseq_dds, normalized=TRUE)\n"
	rscript += "colnames(co) <- colnames(counts)\n"
	rscript += "write.table(factors, file='size_factors.tsv', sep='\\t', quote=F, row.names=F)\n"
	rscript += "write.table(co, file='normalised_counts.tsv', sep='\\t', quote=F)\n"
	return rscript

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Prints DESEQ2 normalised counts and also size factors\n')
	parser.add_argument('-c','--config', help='Config file containing [Conditions], please see documentation for usage!', required=False)
	parser.add_argument('-i','--input', help='Counts matrix',required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])

	#Read design matrix and create list of conditions and directories
	conditions = ConfigSectionMap("Conditions", Config)
	create_design_for_R(conditions)
	
	create_design_for_R(conditions)
	rscript = print_norm_counts(args["input"])
	run_rcode(rscript, "get_counts.R")