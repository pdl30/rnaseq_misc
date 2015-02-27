#!/usr/bin/python

########################################################################
# 12 Jan 2015
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import re, os, sys
import argparse
import tempfile
import pysam
import ConfigParser
from collections import defaultdict

def create_design_for_R(conditions):
	output = open("tmp_design.txt", "w")
	output.write("countfile\tcondition\n"),
	for key in sorted(conditions.keys()):
		output.write("{}\t{}\n".format(key, conditions[key]))
	output.close()

def monocle(conditions):
	rscript = "suppressPackageStartupMessages(library('Biobase'))\n"
	rscript += "suppressPackageStartupMessages(library('knitr'))\n"
	rscript += "suppressPackageStartupMessages(library('reshape2'))\n"
	rscript += "suppressPackageStartupMessages(library('ggplot2'))\n"
	rscript += "suppressPackageStartupMessages(library('monocle'))\n"
	rscript += "pdata <- read.table('tmp_design.txt', header=T)\n"
	rscript += "pd <- new('AnnotatedDataFrame', data = pdata)\n"
	#Need to run DESEQ2 first to get normalised counts and differentially expressed genes for each of the conditions found in the sets
	rscript += "x <- read.table('total_normalised_matrix.tsv')\n"
	rscript += "HSMM <- newCellDataSet(as.matrix(t(wt)), phenoData = wtpd)\n"
	rscript += "HSMM <- detectGenes(HSMM, min_expr = 0.1)\n"
	rscript += "expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 50))\n"
	rscript += "L <- log(exprs(HSMM[expressed_genes,]))\n"
	rscript += "melted_dens_df <- melt(t(scale(t(L))))\n"
	rscript += "pdf("monocle_report.pdf")\n"
	rscript += "qplot(value, geom='density', data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') + xlab('Standardized log(TPM)') + ylab('Density')\n"




def join_counts(idict):
	data = defaultdict(list)
	output = open("sample_spreadsheet.tsv", "w")
	output.write("ID"),
	for count in sorted(idict):
		output.write("\t{}".format(count)),
		with open(count) as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				if word[0].startswith("__"):
					pass
				else:
					data[word[0]].append(word[1])
	output.write("\n"),
	for key2 in sorted(data):
		data2 = data[key2]
		output.write(key2+"\t" + "\t".join(data2) + "\n"),
	output.close()
	return "sample_spreadsheet.tsv" #NEED TO 

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

def main():
	parser = argparse.ArgumentParser(description='Runs monocle, creating spanning tree.\n')
	parser.add_argument('-c', '--config.ini', help='[Conditions] contains the counts files as keys and condition as value', required=True)
	parser.add_argument('-s', '--spread', help='Spreadsheet, names must be same as config.ini, not required', required=False)
	parser.add_argument('-o', '--output', help='Output file')
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])

	conditions = ConfigSectionMap(Config, "Conditions")
	if args["spread"]:
		spread = args["spread"]
	else:
		spread = join_counts(conditions)