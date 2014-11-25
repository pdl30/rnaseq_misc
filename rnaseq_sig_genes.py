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

def get_gfold_sigs(ifile, cutoff):
	col_index = {}
	sig_genes = {}
	data = {}
	with open(ifile) as f:
		header = next(f)
		header = header.rstrip()
		hwords = header.split("\t")
		c = 0
		for i in hwords:
			m = re.search("Gfold LFC", i) #Searches for LFC columns and gets index
			if m:
				col_index[c] = 1
			c += 1
		c2 = 0
		for key in sorted(col_index): #Prints previous counts names
			if c2 == 0:
				start_key = key
			c2 += 1
		print "Ensembl_ID",
		for q in range(1, len(hwords)):
			print "\t{}".format(hwords[q]),
		print "\n",
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			for i, w in enumerate(word):
				if i in col_index:
					if float(w) < 0:
						y = abs(float(w)) #Converts neg to pos and apply cutoff
						if y > float(cutoff):
							data[word[0]] = line
							counts = []
							for q in range(1, start_key):
								counts.append(word[q])
							sig_genes[word[0]] = counts
					else:
						if float(w) > float(cutoff):
							data[word[0]] = line
							counts = []
							for q in range(1, start_key):
								counts.append(word[q])
							sig_genes[word[0]] = counts
	for gene in sorted(data):
		print "{}\n".format(data[gene]),

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Reads RNA-seq table and selects significant genes\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	deseq_parser = subparsers.add_parser('deseq', help="Input is DESEQ2 table. Not working yet")
	gfold_parser = subparsers.add_parser('gfold', help="Input is GFOLD table")
	deseq_parser.add_argument('-i', '--input', help='Input file', required=True)
	deseq_parser.add_argument('-p', '--pval', help='Cutoff applied to both Padj. default=0.05', default=0.05, required=False)

	gfold_parser.add_argument('-i', '--input', help='Input file', required=True)
	gfold_parser.add_argument('-p', '--pval', help='Cutoff applied to both + and - lfc. default=1', default=1, required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	if args["subparser_name"] == "deseq":
		pass
	elif args["subparser_name"] == "gfold":
		get_gfold_sigs(args["input"], args["pval"])
