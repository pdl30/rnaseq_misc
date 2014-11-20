#!/usr/bin/python

########################################################################
# 09 Oct 2014
# Patrick Lombard, Centre for Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import os, re, sys
import argparse
import subprocess

def rename(idir, f_type, prompt, outname=None):
	ifiles = [f for f in os.listdir(idir)]
#	dont_touch =["tophat_report", "logs", "prep_reads", "junctions", "insertions", "deletions"]
	if outname:
		for ifile in ifiles:
			if ifile.endswith(f_type):
				if prompt:
					Join = raw_input('Would you like to move {} to {}.BED?\n'.format(ifile, outname)).lower()
					if Join == "yes" or Join == "y":
						subprocess.call(["mv", ifile, "{}.BED".format(outname)])
					else:
						print "I won't do anything!\n",
				else:
					subprocess.call(["mv", ifile, "{}.BED".format(outname)])
	else:
		name = os.path.basename(os.path.normpath(idir))
		print name
		for ifile in ifiles:
			if ifile.endswith(f_type):
				if prompt:
					Join = raw_input('Would you like to move {} to {}.BED?\n'.format(ifile, name)).lower()
					if Join == "yes" or Join == "y":
						subprocess.call(["mv", ifile, "{}.BED".format(outname)])
					else:
						print "I won't do anything!\n",
				else:
					subprocess.call(["mv", ifile, "{}.BED".format(outname)])

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Takes HTseq output and creates a fake gfold file.\n ')
	parser.add_argument('-i', '--input', help='Input directory containing files processed by pyrna_align.py', required=True)
	parser.add_argument('-t', '--type', help='Type of file to rename, can do BED/bw/bam etc.', required=False)
	parser.add_argument('-p', '--prompt', action="store_true", help='Prompt before changing', required=False)
	parser.add_argument('-o', '--output', help='Prefix for renaming files to, if not provided, will use folder name', required=False)
	args = vars(parser.parse_args())
	rename(args["input"], args["type"], args["prompt"], args["output"])