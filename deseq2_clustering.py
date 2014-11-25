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

def write_deseq(ifile, sample_dict):
	print "==> Running differental expression analysis...\n"
	rscript =  "library(DESeq2); library(gplots); library(RColorBrewer)\n"
	rscript += "pdata <- read.table('tmp_design.txt', header=T)\n"
	#Make sure they match!
	rscript += "counts <- read.table('{}', sep='\\t', header=T, row.names=1)\n".format(ifile)
	rscript += "rnaseq_dds <- DESeqDataSetFromMatrix(countData = counts, colData = data.frame(pdata), design = ~ condition)\n"
	rscript += "rnaseq_dds$condition <- factor(rnaseq_dds$condition, levels=unique(pdata[,3]))\n"
	rscript += "rnaseq_dds <- DESeq(rnaseq_dds)\n"
	rscript += "rld <- rlog(rnaseq_dds)\n"
	rscript += "distsRL <- dist(t(assay(rld)))\n"
	rscript += "mat <- as.matrix(distsRL)\n"
	rscript += "hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)\n"
	rscript += "rownames(mat) <- colnames(mat) <- with(colData(rnaseq_dds), condition)\n"
	rscript += "hc <- hclust(distsRL)\n"
	rscript += "pdf('clustering.pdf'); heatmap.2(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace='none',col = rev(hmcol), margin=c(13, 13))\n"
	rscript += "tmp <- assay(rld); colnames(tmp) <- c('WT_NSCs', 'sox2_null_NSC_1', 'sox2_null_NSC_2'); tmp_len <- ncol(tmp) -1; \n"
	rscript += "for (i in 1:tmp_len) { plot(tmp[,i], tmp[,i+1], xlab=colnames(tmp)[i], ylab=colnames(tmp)[i+1])}\n"
	rscript += "tmp2 <- tmp[which(rowSums(tmp) != 0),]; c <- cor(tmp2, method='spearman'); d <- as.dist(1-c); hr <- hclust(d, method = 'complete', members=NULL)\n"
	rscript += "plotPCA(rld, intgroup=c('condition'))\n"
	rscript += "plot(hr, hang = 0.1);\n"
	rscript += "dev.off()\n"
	return rscript


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Various clustering for RNA-seq experiments using DESEQ2 counts\n')
	parser.add_argument('-c','--config', help='Config file containing parameters, please see documentation for usage!', required=False)
	parser.add_argument('-i','--input', help='Combined counts file from HTSeq or pyrna_count.py',required=True)
	args = vars(parser.parse_args())
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)
	create_design_for_R(conditions)

	rscript = write_deseq(args["input"], conditions) ##Needs changing!!!
	run_rcode(rscript, "deseq2_rcode.R")