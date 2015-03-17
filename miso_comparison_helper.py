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
import itertools
from multiprocessing import Pool

def run_all_possibles(se, cond1, cond2):
#	for e in itertools.combinations(range(1, 5), 2):
	idir1 = "{}_merged_namesort2_miso_{}".format(cond1, se)
	idir2 = "{}_merged_namesort2_miso_{}".format(cond2, se)
	command = "compare_miso --compare-samples {} {} comparisons/".format(idir1, idir2)
	print command
	subprocess.call(command.split())
	command2 = "filter_events --filter comparisons/{0}_vs_{1}/bayes-factors/{0}_vs_{1}.miso_bf --num-inc 1 --num-exc 1 --num-sum-inc-exc 10 --delta-psi 0.20 --bayes-factor 10 --output-dir filtered/".format(idir1, idir2)
	print command2
	subprocess.call(command2.split())

def run_all_possible_fun(args):
	return run_all_possibles(*args)

#def run_sashimi_fun(args):
#	return run_sashimi(*args)

#def run_sashimi(se, cond1, cond2):
#	for e in itertools.combinations(range(1, 5), 2):
#		idir1 = "{}.{}_merged_namesort2_miso_{}".format(cond1, e[0], se)
#		idir2 = "{}.{}_merged_namesort2_miso_{}".format(cond2, e[1], se)
#		command = "sashimi_plot --plot-bf-dist comparisons/{0}_vs_{1}/bayes-factors/{0}_vs_{1}.miso_bf settings.txt --output-dir plots/".format(cond1, cond2)
#		subprocess.call(command.split())


def main():
	parser = argparse.ArgumentParser(description='Runs miso compare-samples and filter_events for all combinations of splicing events provided by miso. Check naming before running\n')
	parser.add_argument('-s', '--cond1', help='Condition1', required=False)
	parser.add_argument('-c', '--cond2', help='Condition2', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	types = ["A3SS", "A5SS", "AFE", "ALE", "MXE", "RI", "SE", "TandemUTR"]

	pool = Pool(8)
	pool.map(run_all_possible_fun, itertools.izip(types, itertools.repeat(args["cond1"]), itertools.repeat(args["cond2"]))) ##Running annotation in parallel
	pool.close()
	pool.join()
	run_all_possibles(args["cond1"], args["cond2"], types)

main()