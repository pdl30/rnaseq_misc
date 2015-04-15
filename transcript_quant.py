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
from collections import defaultdict
import pybedtools 
import ConfigParser
from multiprocessing import Pool
import itertools

def index(fa, index):
	#First download genes and gene predictions ensembl table from UCSC and save as hg19_transcriptome.fa 
	command = "salmon index -t {} -i {}".format(fa, index)
	subprocess.call(command.split())

def quant(bam, fa, index):
	#need to convert bam to fastq
	f = open(os.devnull,"w")
	name = os.path.basename(bam)
	name = name.strip(".bam$")
	command = "bedtools bamtofastq -i {} -fq {}_1.fq -fq2 {}_2.fq".format(bam, name, name)
	subprocess.call(command.split(), stderr=f)
	#command2 = "salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fa -2 reads2.fa -o transcripts_quant"

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

def quant_fun(args):
	return quant(*args)

def main():
	parser = argparse.ArgumentParser(description='Transcript quantification using salmon\n')
	parser.add_argument('-c', '--config', help='Conditions containing BAM files as keys', required=True)
	parser.add_argument('-t','--threads', help='threads, default=1', default=1, required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)

	fa = "/home/patrick/Reference_Genomes/hg19/UCSC/hg19_transcriptome.fa"
	index = "/home/patrick/Reference_Genomes/hg19/UCSC/salmon_index"

	pool = Pool(int(args["threads"]))
	pool.map(quant_fun, itertools.izip(list(conditions.keys()), itertools.repeat(fa), itertools.repeat(index))) ##Running annotation in parallel
	pool.close()
	pool.join()

main()
# -v [ --version ]                      print version string.
 # -h [ --help ]                         produce help message.
 # -l [ --libType ] arg                  Format string describing the library 
 #                                       type.
 # -a [ --alignments ] arg               input alignment (BAM) file(s).
 # -w [ --maxReadOcc ] arg (=200)        Reads "mapping" to more than this many 
 #                                       places won't be considered.
 # -t [ --targets ] arg                  FASTA format file containing target 
 #                                       transcripts.
 # -p [ --threads ] arg (=6)             The number of threads to use 
 #                                       concurrently. The alignment-based 
 #                                       quantification mode of salmon is 
 #                                       usually I/O bound so until there is a 
 #                                       faster multi-threaded SAM/BAM parser to
#                                        feed the quantification threads, one 
#                                        should not expect much of a speed-up 
#                                        beyond ~6 threads.
#  -e [ --useReadCompat ]                [Currently Experimental] : Use the 
#                                        orientation in which fragments were 
#                                        "mapped"  to assign them a probability.
#                                          For example, fragments with an 
#                                        incorrect relative oritenation with 
#                                        respect  to the provided library format
#                                        string will be assigned a 0 
#                                        probability.
#  --incompatPrior arg (=1.0000000000000001e-05)
#                                        This option can only be used in 
#                                        conjunction with --useReadCompat.  It 
#                                        sets the prior probability that an 
#                                        alignment that disagrees with the 
#                                        specified library type (--libType) 
#                                        results from the true fragment origin. 
#                                        Setting this to 0 says that alignments 
#                                        that disagree with the library type 
#                                        should be "impossible", while setting 
#                                        it to 1 says that alignments that 
#                                        disagree with the library type are no 
#                                        less likely than those that do (in this
#                                        case, though, there is no reason to 
                                     #   even use --useReadCompat)
#  --noEffectiveLengthCorrection         Disables effective length correction 
#                                        when computing the probability that a 
#                                        fragment was generated from a 
#                                        transcript.  If this flag is passed in,
                   #                     the fragment length distribution is not#
                    #                    taken into account when computing this #
                     #                   probability.
  #--noFragLengthDist                    [Currently Experimental] : Don't 
                      #                  consider concordance with the learned 
   #                                     fragment length distribution when 
    #                                    trying to determine the probability 
      #                                  that a fragment has originated from a 
     #                                   specified location.  Normally, 
        #                                Fragments with unlikely lengths will be
       #                                 assigned a smaller relative probability
#                                        than those with more likely lengths.  
#                                        When this flag is passed in, the 
#                                        observed fragment length has no effect #
                                     #   on that fragment's a priori 
                                      #  probability.
  #--useErrorModel                       [Currently Experimental] : Learn and 
  #                                      apply an error model for the aligned 
  #                                      reads.  This takes into account the the
  #                                      observed frequency of different types 
  #                                      of mismatches when computing the 
  #                                      likelihood of a given alignment.
  #--numErrorBins arg (=6)               The number of bins into which to divide
  #                                      each read when learning and applying 
 #                                       the error model.  For example, a value 
 #                                       of 10 would mean that effectively, a 
 #                                       separate error model is leared and 
   #                                     applied to each 10th of the read, while
  #                                      a value of 3 would mean that a separate
     #                                   error model is applied to the read 
    #                                    beginning (first third), middle (second
  #                                      third) and end (final third).
 # -f [ --forgettingFactor ] arg (=0.65000000000000002)
   #                                     The forgetting factor used in the 
    #                                    online learning schedule.  A smaller 
     #                                   value results in quicker learning, but 
   #                                     higher variance and may be unstable.  A
      #                                  larger value results in slower learning
       #                                 but may be more stable.  Value should 
        #                                be in the interval (0.5, 1.0].
  #--mappingCacheMemoryLimit arg (=2000000)
         #  #                             If the file contained fewer than this 
          #                              many mapped reads, then just keep the 
         #                               data in memory for subsequent rounds of
        #                                inference. Obviously, this value should
       #                                 not be too large if you wish to keep a 
      #                                  low memory usage, but setting it large 
     #                                   enough to accommodate all of the mapped
    #                                    read can substantially speed up 
   #                                     inference on "small" files that contain
  #                                      only a few million reads.
  #-o [ --output ] arg                   Output quantification directory.
  #-s [ --sampleOut ]                    Write a "postSample.bam" file in the 
        #                                output directory that will sample the 
       #                                 input alignments according to the 
      #                                  estimated transcript abundances. If 
     #                                   you're going to perform downstream 
    #                                    analysis of the alignments with tools 
   #                                     which don't, themselves, take fragment 
 #                                       assignment ambiguity into account, you 
  #                                      should use this output.
 # -u [ --sampleUnaligned ]              In addition to sampling the aligned 
  #                                      reads, also write the un-aligned reads 
  #                                      to "posSample.bam".
  #--biasCorrect                         [Experimental]: Output both 
  #                                      bias-corrected and non-bias-corrected 
  #                                      qunatification estimates.
  #-n [ --numRequiredObs ] arg (=50000000)
     #                                   The minimum number of observations 
    ##                                    (mapped reads) that must be observed 
    #                                    before the inference procedure will 
    ##                                    terminate.  If fewer mapped reads exist
   #                                     in the input file, then it will be read
  #                                      through multiple times.
  #-g [ --geneMap ] arg                  File containing a mapping of 
                         #               transcripts to genes.  If this file is 
                        #                provided Sailfish will output both 
                       #                 quant.sf and quant.genes.sf files, 
                      #                  where the latter contains aggregated 
                     #                   gene-level abundance estimates.  The 
                   ##                     transcript to gene mapping should be 
                  #                      provided as either a GTF file, or a in 
                 #                       a simple tab-delimited format where 
                #                        each line contains the name of a 
               #                         transcript and the gene to which it 
              #                          belongs separated by a tab.  The 
             #                           extension of the file is used to 
            #                            determine how the file should be 
           #                             parsed.  Files ending in '.gtf' or 
          #                              '.gff' are assumed to be in GTF format;
         #                               files with any other extension are 
        #                                assumed to be in the simple format.
