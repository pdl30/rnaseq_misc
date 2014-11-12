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

def htseq_to_gfold(ilist, gfold_file):
	gfold_data = {}
	with open(gfold_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			gfold_data[word[0]] = (word[1], word[3])
	for ifile in ilist:
		name = ifile.strip(".count")
		output = open(name+"_htgold.count", "w")
		with open(ifile) as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				m = re.search("__", word[0])
				if m:
					pass
				else:
					output.write("{}\t{}\t{}\t{}\t{}\n".format(word[0], gfold_data[word[0]][0], word[1], gfold_data[word[0]][1], 10)),
		output.close()

def table_to_gfold(ifile, gfold_file):
	gfold_data = {}
	with open(gfold_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			gfold_data[word[0]] = (word[1], word[3]) #This is gene id, name and gene length

	with open(ifile) as f:
		data = {}
		header = next(f)
		header = header.rstrip()
		head = header.split("\t")

		for i in range(1, len(head)):
			data[head[i]] = {}
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			m = re.search("__", word[0])
			if m:
				pass
			else:
				for i in range(1, len(word)):
					data[head[i]][word[0]] = word[i]
	for key in sorted(data):
		output = open("{}_htgold.tsv".format(key), "w")
		for gene in sorted(data[key]):
			output.write("{}\t{}\t{}\t{}\t{}\n".format(gene, gfold_data[gene][0], data[key][gene], gfold_data[gene][1], 10)),
		output.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Takes HTseq output and creates a fake gfold file.\n ')
	parser.add_argument('-i', '--input', help='Input directory containing files with .count ending', required=False)
	parser.add_argument('-m', '--matrix', help='Alternatively use input matrix. Will output _htgold.tsv files', required=False)
	parser.add_argument('-e', '--ex', help='Example gfold count file processed using same GTF', required=True)
	args = vars(parser.parse_args())
	if args["input"]:
		ifiles = [f for f in os.listdir(args["input"]) if f.endswith(".count")]
		htseq_to_gfold(ifiles, args["ex"])
	elif args["matrix"]:
		table_to_gfold(args["matrix"], args["ex"])
