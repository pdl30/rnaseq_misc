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
import ConfigParser

#No point, way to variable!! 
def create_staging_table(accession, out):
	output = open(out, "w")
	header = """id\tCT_General\tCT_Subtype\tbto\tfactor\tfactorGeneId\tdetails\tGSE\tGSM\tgroupname\trepository\texp_type\tfilename\tcontrol\tSRX\tpvaluemerged\t
		paired\tmapping_genome\ttotal_raw_reads\ttrimmed\treads_after_trimming\tuniquely_mappable_reads\tFinal_peak_number\taligner_used\traw_data_code\tIn_house_code\t
		Completion_date\tAdded_to_compendium\tstatus\tprocessing_comments\talignment_report\torganism\tsubmitter\n"""
	output.write(header),
	command = "wget -c -nv -q -O tmp.txt http://www.ebi.ac.uk/arrayexpress/files/{0}/{0}.sdrf.txt".format(accession)
	command = "wget -c -nv -q -O tmp2.txt http://www.ebi.ac.uk/arrayexpress/files/{0}/{0}.idf.txt".format(accession)
	subprocess.call(command, shell=True)

	f = open("tmp.txt", "r")
	lines = f.readlines()
	header = lines[0].rstrip()
	head = header.split("\t")
	cell_type =  [i for i, x in enumerate(head) if x.startswith("Characteristics") and "cell type" in x]
	cell_line =  [i for i, x in enumerate(head) if x.startswith("Characteristics") and "cell line" in x]
	organism =  [i for i, x in enumerate(head) if x.startswith("Characteristics") and "organism" in x]
	library_layout = [i for i, x in enumerate(head) if x.startswith("Comment") and "LIBRARY_LAYOUT" in x]
	library_type = [i for i, x in enumerate(head) if x.startswith("Comment") and "LIBRARY_STRATEGY" in x]
	names = [i for i, x in enumerate(head) if "Assay Name" in x]

	for i in xrange(1, len(lines)):
		line = lines[i]
		line = line.rstrip()
		word = line.split("\t")
		#print cell_type, cell_line, organism, library_layout, library_type
		print word[cell_type[0]], word[cell_line[0]], word[organism[0]], word[library_layout[0]], word[library_type[0]]

	with open("tmp2.txt") as f:
		for line in f:
			if line.startswith("Protocol Description"):
				description = re.sub("Protocol Description\t", "", line)


def download_ebi(accession):
	#Example EBI Format
	#Source Name	Comment[ENA_SAMPLE]	Material Type	Provider	Characteristics[organism]	Characteristics[specimen with known storage state]	
	#Characteristics[cell line]	Characteristics[strain]	Characteristics[cell type]	
	#Protocol REF	Protocol REF	Protocol REF	Protocol REF	Extract Name	Material Type	
	#Comment[LIBRARY_LAYOUT]	Comment[LIBRARY_SOURCE]	Comment[LIBRARY_STRATEGY]
	command = "wget -c -nv -q -O tmp.txt http://www.ebi.ac.uk/arrayexpress/files/{0}/{0}.sdrf.txt".format(accession)
	subprocess.call(command, shell=True)
	f = open("tmp.txt", "r")
	lines = f.readlines()
	header = lines[0].rstrip()
	head = header.split("\t")
	names = [i for i, x in enumerate(head) if x.startswith("Comment") and "ENA_SAMPLE" in x]
	links = [i for i, x in enumerate(head) if x.startswith("Comment") and "FASTQ_URI" in x]
	old_path = os.getcwd()
	for i in xrange(1, len(lines)):
		line = lines[i]
		line = line.rstrip()
		word = line.split("\t")
		if not os.path.isdir(word[names[0]]):
			os.mkdir(word[names[0]])
		os.chdir(word[names[0]])
		link = word[links[0]]
		fastq = os.path.basename(link)
		command = "wget -c -nv -q {0}".format(link)
		subprocess.call(command.split())
		command2 = "gunzip {0}".format(fastq)
		subprocess.call(command, shell=True)
		os.chdir(old_path)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Downloads samples from EBI\n')

	parser.add_argument('-a', '--accession', help='Arrayexpress accession number', required=True )
	#parser.add_argument('-s', action='store_true', help='This will create a staging table and won\'t process the samples', required=False )
	#parser.add_argument('-o', '--output', help='Output staging table', required=False )
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	download_ebi(args["accession"])
