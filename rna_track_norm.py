#!/usr/bin/python

########################################################################
# 15 May 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, os, re
import pybedtools
import pysam
import argparse
import operator
import pkg_resources
import pychiptools

##Must include scaling!
def genomeCoverage(name, house):
	print "==> Converting bed to bedGraph...\n"
	inbed = pybedtools.BedTool(name+"_ucsc.BED")
	outcov = inbed.genome_coverage(bg=True, genome='mm10', scale=house)
	output = name+"_house.bedGraph"
	outcov.saveas(output)
	return output

def bedgraphtobigwig(bedgraph, chrom):
	bw = re.sub(".bedGraph$", ".bw", bedgraph)
	print "==> Converting bedGraph to bigWig...\n"
	command = ["bedGraphToBigWig", bedgraph, chrom, bw]
	subprocess.call(command)

def normalise_to_housekeeper(count_file):
	#print "==> Normalising to Housekeeper...\n"
	with open(count_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if word[0] == "ENSMUSG00000057666": #Gapdh, substitute with what you want to use. REmove from production?
				housekeeper = int(word[1])	
	return housekeeper

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Processes RNA-seq samples to bigWig tracks.\nIf Tophat2 is specified, this will pull out the uniquely mapped reads\nOtherwiseit is assumed that the bam file is already uniquely mapped!")
	parser.add_argument('-i', '--input', help='BED file in UCSC format', required=True)
	parser.add_argument('-g', '--genome', help='Genome the samples are aligned to, options include mm10/mm9/hg19', required=True)
	parser.add_argument('-a', '--house', help='Housekeeper normalisation. Input file is HTSEQ-count file containing gene for normalisation on first line', required=False)
	args = vars(parser.parse_args())

	chrom = pkg_resources.resource_filename('pyrnatools', 'data/{}.chrom.sizes'.format(args["genome"]))
	if not os.path.isfile(chrom):
		raise Exception("Unsupported Genome!")
	
	name = re.sub("_ucsc.BED$", "", args["input"])

	house = normalise_to_housekeeper(args["house"])
	scale = float(1000)/int(house) #Works and checked
	bedgraph = genomeCoverage(name, house=scale)

	bedgraphtobigwig(bedgraph, chrom)
