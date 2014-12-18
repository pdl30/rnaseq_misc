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
def genomeCoverage(name, genome, house=None, deseq=None, split=False):
	print "==> Converting bed to bedGraph...\n"
	inbed = pybedtools.BedTool(name+"_ucsc.BED")
	if split:
		if house:
			outcov1 = inbed.genome_coverage(bg=True, strand="+", genome=genome, scale=house)
			output1 = name+"_pos_house.bedGraph"
			outcov2 = inbed.genome_coverage(bg=True, strand="-", genome=genome, scale=house)
			output2 = name+"_neg_house.bedGraph"
			output = [output1, output2]
			outcov1.saveas(output1)
			outcov2.saveas(output2)
		elif deseq:
			outcov1 = inbed.genome_coverage(bg=True, strand="+", genome=genome, scale=deseq)
                        output1 = name+"_pos_deseq.bedGraph"
                        outcov2 = inbed.genome_coverage(bg=True, strand="-", genome=genome, scale=deseq)
                        output2 = name+"_neg_deseq.bedGraph"
                        output = [output1, output2]
			outcov1.saveas(output1)
                        outcov2.saveas(output2)
	else:
		if house:
			outcov = inbed.genome_coverage(bg=True, genome=genome, scale=house)
			output = name+"_house.bedGraph"
		elif deseq:
			outcov = inbed.genome_coverage(bg=True, genome=genome, scale=deseq)
			output = name+"_deseq2.bedGraph"
		outcov.saveas(output)
	return output

def bedgraphtobigwig(bedgraph, chrom, split=False):
	if split:
		for bedg in bedgraph:
			bw = re.sub(".bedGraph$", ".bw", bedg)
		        print "==> Converting bedGraph to bigWig...\n"
       			command = ["bedGraphToBigWig", bedg, chrom, bw]
			subprocess.call(command)
	else:
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
	parser = argparse.ArgumentParser(description="Processes RNA-seq samples to bigWig tracks")
	parser.add_argument('-i', '--input', help='BED file in UCSC format', required=True)
	parser.add_argument('-g', '--genome', help='Genome the samples are aligned to, options include mm10/mm9/hg19', required=True)
	parser.add_argument('-a', '--house', help='Housekeeper normalisation. Input file is HTSEQ-count file containing gene for normalisation on first line', required=False)
	parser.add_argument('-d', '--deseq2', help='DESEQ2 sizeFactor normalisation')	
	parser.add_argument('-s', '--split', help='Splits the bigwig by strand', action='store_true')
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	chrom = pkg_resources.resource_filename('pyrnatools', 'data/{}.chrom.sizes'.format(args["genome"]))
	if not os.path.isfile(chrom):
		raise Exception("Unsupported Genome!")
	
	name = re.sub("_ucsc.BED$", "", args["input"])
	name = re.sub(".BED$", "", name)
	if args["house"]:
		house = normalise_to_housekeeper(args["house"])
		scale = float(1000)/int(house) #Works and checked
		bedgraph = genomeCoverage(name, args["genome"], house=scale, split=args["split"])
	elif args["deseq2"]:
		sizeF = 1/float(args["deseq2"])
		bedgraph = genomeCoverage(name, args["genome"], deseq=sizeF, split=args["split"])
	bedgraphtobigwig(bedgraph, chrom, args["split"])
