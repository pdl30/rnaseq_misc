#!/usr/bin/python

########################################################################
# 12 Jan 2015
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, re, os
import argparse
import ConfigParser
from multiprocessing import Pool
import itertools

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

def seqgsea_count(sample, path, gtf, paired, orientation):
	bam_name = os.path.basename(sample)
	output = re.sub(".bam$", "_seqgsea.count", bam_name)
	if paired:
		p = "yes"
	else:
		p = "no"
	command = "python {} -b yes -p {} -s {} {} {} {}".format(path, p, orientation, gtf, sample, output)
	subprocess.call(command.split())

def run_seqgsea(conditions, comp1, comp2):
	#Current problem is that this takes fecking ages
	rscript = 'library("SeqGSEA")\n'
	rscript += "pdata <- read.table('tmp_design.txt', header=T)\n"
	rscript += "counts1 <- pdata[which(pdata[,2] == '{}'),1]\n".format(comp1)
	rscript += "counts2 <- pdata[which(pdata[,2] == '{}'),1]\n".format(comp2)
	rscript += "RCS <- loadExonCountData(as.character(counts1), as.character(counts2))\n"
	rscript += "RCS <- exonTestability(RCS, cutoff=5)\n"
	rscript += "geneTestable <- geneTestability(RCS)\n"
	rscript += "RCS <- subsetByGenes(RCS, unique(geneID(RCS))[ geneTestable ])\n"
	rscript += "geneIDs <- unique(geneID(RCS))\n"
	rscript += "RCS <- estiExonNBstat(RCS)\n"
	rscript += "RCS <- estiGeneNBstat(RCS)\n"
	rscript += "perm.times <- 1000\n"
	rscript += "permuteMat <- genpermuteMat(RCS, times=perm.times)\n"
	rscript += "RCS <- DSpermute4GSEA(RCS, permuteMat)\n"

def seqgsea_count_fun(args):
	return seqgsea_count(*args)

def spliceR(idir, gtf, genome):
	#Uses cuffdiff directories. Not working currently, giving errors
	rscript = "library(spliceR)\n"
	rscript += "cuffDB <- readCufflinks(dir={},gtf={},genome='{}')\n"
	rscript += "cuffDB_spliceR <- prepareCuff(cuffDB)\n"
	rscript += "myTranscripts <- transcripts(cuffDB_spliceR); myExons <- exons(cuffDB_spliceR); conditions(cuffDB_spliceR)\n"
	#rscript += "cuffDB_spliceR_filtered <- preSpliceRFilter(cuffDB_spliceR,filters=c('expressedIso', 'isoOK', 'expressedGenes', 'geneOK'))\n"
	rscript += "mySpliceRList <- spliceR(cuffDB_spliceR, compareTo='preTranscript', filters=c('expressedGenes','geneOK', 'isoOK', 'expressedIso', 'isoClass'))\n"
	rscript += "ucscCDS <- getCDS(selectedGenome='hg19', repoName='UCSC'); require('BSgenome.Hsapiens.UCSC.hg19', character.only = TRUE)\n"
	rscript += "PTCSpliceRList <- annotatePTC(cuffDB_spliceR, cds=ucscCDS, Hsapiens, PTCDistance=50)\n"
	rscript += "generateGTF(mySpliceRList, filters=c('geneOK', 'isoOK', 'expressedGenes', 'expressedIso'), scoreMethod='local', useProgressBar=F)\n"
	rscript += "mySpliceRList <- spliceRPlot(mySpliceRList, evaluate='nr_transcript_pr_gene')\n"
	rscript += "mySpliceRList <- spliceRPlot(mySpliceRList, evaluate='mean_AS_transcript', asType='ESI')\n"
	return rscript

def create_design_for_R(idict):
	output = open("tmp_design.txt", "w")
	output.write("sampleName\tfileName\tcondition\n"),
	for key in sorted(idict.keys()):
		bam_name = os.path.basename(sample)
		name = re.sub(".bam$", "", bam_name)
		count = re.sub(".bam$", "_dexseq.count", bam_name)
		output.write("{}\t{}\t{}\n".format(name, count, idict[key]))
	output.close()

def reverse_dict(idict):
	inv_map = {}
	for k, v in idict.iteritems():
		inv_map[v] = inv_map.get(v, [])
		inv_map[v].append(k)
	return inv_map

def main():
	parser = argparse.ArgumentParser(description='Overview of a few programs for Splicing analysis\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")

	splice_parser = subparsers.add_parser('spliceR', help="Runs spliceR")
	splice_parser.add_argument('-c','--config', help='Config file containing cuffdiff directories as keys', required=True)
	splice_parser.add_argument('-g','--gtf', help='GTF file formatted by spliceR', required=True)
	splice_parser.add_argument('-o','--output', help='Output directory', required=True)

	seq_parser = subparsers.add_parser('seqGSEA', help="Runs seqGSEA")
	seq_parser.add_argument('-c','--config', help='Config file containing bam files, please see documentation for usage!', required=True)
	seq_parser.add_argument('-g','--gtf', help='GTF file formatted by prepare_exon_annotation_ensembl.py script', required=True)
	seq_parser.add_argument('-t','--threads', help='threads, default=1', default=1, required=False)
	seq_parser.add_argument('-p', action='store_true', help='Use if samples are paired end. Will find sd and insert size for bam files', required=False)
	seq_parser.add_argument('-o','--orientation', help='Options are yes, no or reverse. Test First!!!', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)

	if args["subparser_name"] == "spliceR": #Not even close to finished
		for key in conditions:
			spliceR(key, args["gtf"], args["output"])

	elif args["subparser_name"] == "seqGSEA":
		path = "/raid/home/patrick/R/x86_64-pc-linux-gnu-library/3.1/SeqGSEA/extscripts"
		count_program = path + "/count_in_exons.py"
		pool = Pool(int(args["threads"]))
		pool.map(seqgsea_count_fun, itertools.izip(list(conditions.keys()), itertools.repeat(count_program), itertools.repeat(args["gtf"]), itertools.repeat(args["p"]),
			itertools.repeat(args["orientation"]))) ##Running annotation in parallel
		pool.close()
		pool.join()
		comparisons = ConfigSectionMap("Comparisons", Config)
		for comp in comparisons:
			c = comparisons[comp].split(",")
			comps = [x.strip(' ') for x in c]
			rscipt = run_seqgsea(conditions, comps[0], comps[1])
			run_rcode(rscript, "dexseq.R")
			
main()
