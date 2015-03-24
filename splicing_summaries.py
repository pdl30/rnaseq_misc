	#!/usr/bin/python

########################################################################
# 9 Jan 2015
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import argparse
import os, re, sys
import subprocess
import HTSeq
from collections import defaultdict

def read_dexseq():
	dif1 = defaultdict(list)
	dif2 = defaultdict(list)
	with open("dexseq/75K.DIF_vs_CTR.DIF_dexseq.tsv") as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if len(word) > 5:
				if word[5] == "NA":
					pass
				else:
					pval = float(word[5])
					if pval <= 0.05:
						exon = word[0].split(":")
						dif1[word[1]].append("{}:{}".format(pval, exon[1]))
	with open("dexseq/75K.UND_vs_CTR.UND_dexseq.tsv") as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if len(word) > 5:
				if word[5] == "NA":
					pass
				else:
					pval = float(word[5])
					if pval <= 0.05:
						exon = word[0].split(":")
						dif2[word[1]].append("{}:{}".format(pval, exon[1]))
	return dif1, dif2

def read_mats():
	dif1 = {}
	dif2 = {}
	events = ["A3SS", "A5SS", "RI", "MXE", "SE"]

	for event in events:
		with open(os.path.join("mats/75K.DIF_vs_CTR.DIF/MATS_output/", "{}.MATS.ReadsOnTargetAndJunctionCounts.txt".format(event))) as f:
			head = next(f)
			header = head.rstrip().split("\t")
			c = 0
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				try:
					if event == "MXE":
						pval = float(word[20])
					else:
						pval = float(word[18])
				except ValueError:
					print c
				c += 1
				name = word[1].strip("\"")
				dif1[name] = [pval, event]

		with open(os.path.join("mats/75K.UND_vs_CTR.UND/MATS_output/", "{}.MATS.ReadsOnTargetAndJunctionCounts.txt".format(event))) as f:
			next(f)
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				try:
					if event == "MXE":
						pval = float(word[20])
					else:
						pval = float(word[18])
				except ValueError:
					print c
				name = word[1].strip("\"")
				dif2[name] = [pval, event]
	return dif1, dif2

def read_cuff(conv):
	dif1 = {}
	dif2 = {}
	with open("cuff/75K.DIF_vs_CTR.DIF/splicing.diff") as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			#if word[13] == "yes":
			pval = float(word[11])
			if word[2] in conv:
			#print word[2], conv[word[2]]
				dif1[conv[word[2]]] = pval
	with open("cuff/75K.UND_vs_CTR.UND/splicing.diff") as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			#if word[13] == "yes":
			pval = float(word[11])
			if word[2] in conv:
				dif2[conv[word[2]]] = pval
	return dif1, dif2

def read_anno():
	conv = {}
	inv_conv = {}
	gtffile = HTSeq.GFF_Reader( "../annotation/Homo_sapiens.GRCh37.74.ucsc.gtf" )
	for feature in gtffile:
		conv[feature.attr["gene_id"]] = feature.attr["gene_name"]
		inv_conv[feature.attr["gene_name"]] = feature.attr["gene_id"]
	return conv, inv_conv

def write_results(dex, mats, cuf, miso, conv, out):
	output = open(out, "w")
	output.write("Ensembl_ID\tGene_Name\tDEXSEQ pval\tMATS pval\tMATS event\tCUFFDIFF pval\tMISO bayes\tMISO event\tMISO significant Y/N\n")
	for key in conv:
		name1 = dex.get(key, None)
		name2 = mats.get(key, ["NA", "NA"])
		name3 = cuf.get(key, "NA")
		name4 = miso.get(key, ["NA", "NA", "NA"])
		if name1:
			output.write("{}\t{}\t".format(key, conv[key])),
			output.write(",".join(name1)),
		#	print name4
			output.write("\t{}\t{}\t{}\t{}\t{}\t{}\n".format(name2[0], name2[1], name3, name4[0], name4[1], name4[2])),
		else:
			output.write("{}\t{}\tNA\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key, conv[key], name2[0], name2[1], name3, name4[0], name4[1], name4[2])),

def miso_conversion():
	events = ["A3SS", "A5SS", "RI", "MXE", "SE"]
	conv = {}
	for eve in events:
		conv[eve] = {}
		with open("miso/hg19/{}.hg19.gff3".format(eve)) as f:
			next(f)
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				if word[2] == "gene":
					info = word[8].split(";")
					name=  info[0].strip("^Name=")
					ens = info[3].strip("^ensg_id=")
					conv[eve][name] = ens
	return conv

def read_miso(id_conv):
	dif1 = {}
	dif2 = {}
	events = ["A3SS", "A5SS", "RI", "MXE", "SE"]
	for eve in events:
		dif1[eve] = {}
		name = "75K.DIF_merged_namesort2_miso_{0}_vs_CTR.DIF_merged_namesort2_miso_{0}".format(eve)
		with open("miso/comparisons/{0}/bayes-factors/{0}.miso_bf".format(name)) as f:
			next(f)
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				bayes = float(word[8])
				if word[0] in id_conv[eve]:
					name = id_conv[eve][word[0]]
					if abs(float(word[7])) >= 0.2 and bayes >= 10:
						dif1[name] = [bayes, eve, "Yes"]
					else:
						dif1[name] = [bayes, eve, "No"]
		name = "75K.UND_merged_namesort2_miso_{0}_vs_CTR.UND_merged_namesort2_miso_{0}".format(eve)
		with open("miso/comparisons/{0}/bayes-factors/{0}.miso_bf".format(name)) as f:
			next(f)
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				bayes = float(word[8])

				if word[0] in id_conv[eve]:
					name = id_conv[eve][word[0]]
					if abs(float(word[7])) >= 0.2 and bayes >= 10:
						dif2[name] = [bayes, eve, "Yes"]
					else:
						dif2[name] = [bayes, eve, "No"]
	return dif1, dif2

def main():
#	parser = argparse.ArgumentParser(description='Create summary tables from the different splicing results\n')
#	parser.add_argument('-c', '--config', help='Conditions containing GSM keys', required=False)
#	if len(sys.argv)==1:
#		parser.print_help()
#		sys.exit(1)
#	args = vars(parser.parse_args())


	anno_conv, inv_conv = read_anno()
	dex1, dex2 = read_dexseq() #Both these give ENSG ids
	mats1, mats2 = read_mats()
	cuf1, cuf2 = read_cuff(inv_conv) #This gives gene names, infuriating!!
	miso_conv = miso_conversion()
	miso1, miso2 = read_miso(miso_conv)
	write_results(dex1, mats1, cuf1, miso1, anno_conv, "75K.DIF_vs_CTR.DIF_summary.tsv")
	write_results(dex2, mats2, cuf2, miso2, anno_conv, "75K.UND_vs_CTR.UND_summary.tsv")
main()