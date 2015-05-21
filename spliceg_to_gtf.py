#!/usr/bin/python

########################################################################
# 15 May 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import re, os, sys
import argparse


def read_model(model):
	data = {}
	with open(model) as f:
		for line in f:
			iso = None
			ID = None
			Parent = None
			line = line.rstrip()
			word = line.split("\t")
			if len(word) < 8:
				pass
			else:
				info = word[8].split(";")
				for i in info:
					d = i.split("=")
					if d[0] == "ID":
						ID = d[1]
					elif d[0] == "Parent":
						Parent = d[1]
					elif d[0] == "Isoforms":
						iso = d[1]
				if word[2] == "parent":
					data[iso] = {}
					data[iso]["chr"] = word[0]
					data[iso]["start"] = [word[3]]
					data[iso]["end"] = [word[4]]
					data[iso]["strand"] = word[6]
				elif word[2] == "child":
					isoforms = iso.split(",")
					for i in isoforms:
						if i in data:
							data[i]["start"].append(word[3])
							data[i]["end"].append(word[4])
	return data

def write_gtf(data, gene, output):	
	for transcript in data.keys():
		if data[transcript]["strand"] == "+":
			starts = [int(x) for x in data[transcript]["start"]]
			starts.sort(key=int)
			ends = [int(x) for x in data[transcript]["end"]]
			ends.sort(key=int)
			for i, s in enumerate(starts):
				j = i + 1
				info_string = '''exon_number "{}"; gene_id "{}"; gene_name "{}"; transcript_id "{}"'''.format(j, gene, gene, transcript)
				output.write("{}\tSpliceGrapher\texon\t{}\t{}\t.\t+\t.\t{}\n".format(data[transcript]["chr"], s, ends[i], info_string))
		elif data[transcript]["strand"] == "-":
			starts = [int(x) for x in data[transcript]["start"]]
			starts.sort(key=int)
			ends = [int(x) for x in data[transcript]["end"]]
			ends.sort(key=int)
			for i, s in enumerate(starts):
				j = i + 1
				info_string = '''exon_number "{}"; gene_id "{}"; gene_name "{}"; transcript_id "{}"'''.format(j, gene, gene, transcript)
				output.write("{}\tSpliceGrapher\texon\t{}\t{}\t.\t-\t.\t{}\n".format(data[transcript]["chr"], s, ends[i], info_string))

def main():
	parser = argparse.ArgumentParser(description='Converts gene model to GTF format\n')
	parser.add_argument('-i', '--input', help='Input list of gene model forms', required=True)
	parser.add_argument('-o', '--output', help='GTF output', required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	output = open(args["output"], "w")
	with open(args["input"]) as f:
		for line in f:
			line = line.rstrip()
			name = os.path.basename(line)
			path_to_file = "forms/{}".format(name)
			if os.path.isfile(path_to_file):
				gene = re.sub(".gff", "", name)
				data = read_model(path_to_file)
				write_gtf(data, gene, output)
			else:
				print path_to_file, "is missing!"
	output.close()
main()