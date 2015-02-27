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
import HTSeq
import itertools
import argparse
from multiprocessing import Process, Queue
import pysam

def convert_bam_bed(bam_file, paired, que):
	count = 0
	filtered_bam = pysam.view( "-bq 50", bam_file) ##Filters for uniquely aligned reads!
	for read in filtered_bam:
		count += 1 
	if paired:
		count /= 2
	que.put(count)

def annotate_sam(bam_file, gtf_file, stranded, outfile=None):
	print "==> Counting sam file...\n"
	if outfile:
		htout = open(outfile, "w")
	else:
		count_file = re.sub(".bam$", ".count", bam_file)
		htout = open(count_file,"w")
	command = ["htseq-count", "--mode=union", "--stranded={}".format(stranded), "--quiet", "-f", "bam", bam_file,  gtf_file]
	subprocess.call(command, stdout=htout)
	return count_file

def process_gtf(igtf, que):
	gtf_info = {}
	gtffile = HTSeq.GFF_Reader( igtf )
	for feature in gtffile:
		if feature.type == "exon":
			size = feature.iv.end - feature.iv.start
			gtf_info[feature.name] = size
	que.put(gtf_info)

def calculate_fpkm(gtf_info, count_file, library_size, output_file):
	counts = {}
	with open(count_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if word[0].startswith("__"):
				pass
			else:
				counts[word[0]] = word[1]
	output = open(output_file, "w")
	for c in counts:
		fp = (1000000000 * counts[c])
		km = (library_size * gtf_info[c])
		#print fp, km
		fpkm = float(fp)/float(km)
		output.write("{}\t{}\n".format(c, fpkm)),

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Calculates FPKM from HTSEQ-counts\n')
	parser.add_argument('-i','--input', help='Input BAM file',required=True)
	parser.add_argument('-p', action='store_true', help='Use if sample is paired end', required=False)
	parser.add_argument('-c', '--count', help='HTSEQ-count file', required=False)
	parser.add_argument('-g','--gtf', help='GTF file',required=True)
	parser.add_argument('-s','--stranded', help='Option for HTSeq, default=yes', default="yes", required=False)
	parser.add_argument('-o','--output', help='Output file',required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
#	if args["input"]:
#		count_file = annotate_sam(args["input"], args["gtf"], args["stranded"], None)
#	elif args["count"]:
#		count_file = args["count"]

	queue1 = Queue() #create a queue object
	queue2 = Queue() #create a queue object
	p = Process(target= process_gtf, args= (args["gtf"],queue1))
	p2 = Process(target= convert_bam_bed, args=(args["input"], args["p"], queue2))
	p.start()
	p2.start()
	gtf_pos = queue1.get() #and we're getting return value: 20
	size = queue2.get()
	p.join()
	p2.join()
	#size = convert_bam_bed
	calculate_fpkm(gtf_pos, args["count"], size, args["output"])