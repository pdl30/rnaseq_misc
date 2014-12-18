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
import argparse
import pysam


def convert_bam_bed(name, paired):
	count = 0
	filtered_bam = pysam.view( "-bq 50", name+".bam") ##Filters for uniquely aligned reads!
	for read in filtered_bam:
		count += 1 
	if paired:
		count /= 2
	return count


def main():
	parser = argparse.ArgumentParser(description="Prints library size from bam file")
	parser.add_argument('-i', '--input', help='Bam file from aligner etc.', required=False)
	parser.add_argument('-p', action='store_true', help='Use if samples are paired end.', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	name = re.sub(".bam$", "", args["input"])
	unique_reads = convert_bam_bed(name, args["p"])
	print "{}\t{}\n".format(name, unique_reads),
main()