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
		for ifile in sorted(ifiles):
			if ifile.endswith(f_type):
				if prompt:
					Join = raw_input('Would you like to move {} to {}{}?\n'.format(ifile, outname, f_type)).lower()
					if Join == "yes" or Join == "y":
						subprocess.call(["mv", idir+"/"+ifile, "{}/{}{}".format(idir, outname, f_type)])
					else:
						print "I won't do anything!\n",
				else:
					subprocess.call(["mv", idir+"/"+ifile, "{}/{}{}".format(idir, outname, f_type)])
	else:
		name = os.path.basename(os.path.normpath(idir))
		for ifile in sorted(ifiles):
			if ifile.endswith(f_type):
				if prompt:
					Join = raw_input('Would you like to move {} to {}{}?\n'.format(ifile, name, f_type)).lower()
					if Join == "yes" or Join == "y":
						subprocess.call(["mv", idir+"/"+ifile, "{}/{}{}".format(idir, name, f_type)])
					else:
						print "I won't do anything!\n",
				else:
					subprocess.call(["mv", idir+"/"+ifile, "{}/{}{}".format(idir, name, f_type)])

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Takes HTseq output and creates a fake gfold file.\n ')
	parser.add_argument('-i', '--input', help='Input directory containing files processed by pyrna_align.py', required=True)
	parser.add_argument('-t', '--type', help='Type of file to rename, can do .BED/.bw/.bam etc.', required=False)
	parser.add_argument('-p', '--prompt', action="store_true", help='Prompt before changing', required=False)
	parser.add_argument('-o', '--output', help='Prefix for renaming files to, if not provided, will use folder name', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	rename(args["input"], args["type"], args["prompt"], args["output"])