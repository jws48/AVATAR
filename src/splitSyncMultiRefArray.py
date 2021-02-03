#!/usr/bin/env python

import os
import subprocess
import argparse
import sys

"""
This script is called as a single sample submission or Slurm array job, and takes in
an integer for the number of reference targets, a comma separated list of the reference 
fasta files, path to R1 reads, path to R2 reads, the fasta files containing the forward 
primers to be trimmed, the fasta file containing the reverse primers to be trimmed, 
and the output directory. 

Synchronized fastq paired-end read files are output into the specified directory.

*************************************************************************************
Forward and reverse reads must be stored in separate directories.
*************************************************************************************

command: ./splitSyncMultiRefArray.py -numTargets <numTargs> -refSeqs <ref_seqs> -r1 <read1> -r2 <read2> -fPrimers <for> -rPrimers <rev> -o <out_dir>
"""

#Functions--------------------------------------------------------------------------
def getRefs(numRefs, refs):
	i = 0
	refList = refs.split(',')
	for i in range(refList):
		print(refList[i])


#Main-------------------------------------------------------------------------------
if __name__ == "__main__":

	description = "Splits and synchronizes paired-end amplicon sequence data by target" \

	parser = argparse.ArgumentParser(prog="Variant caller", description=description)
	# user-defined parameters
	parser.add_argument("-numTargets", "--numTargs",
		type=int,
		default=None,
		required=True,
		help="Number of amplicon targets")

	parser.add_argument("-refSeqs", "--ref_seqs",
		type=str,
		default=None,
		required=True,
		help="Comma separated list of paths to reference fasta files")

	parser.add_argument("-r1", "--read1",
		type=str,
		default=None,
		required=True,
		help="Path to forward reads")

	parser.add_argument("-r2", "--read2",
		type=str,
		default=None,
		required=True,
		help="Path to reverse reads")

	parser.add_argument("-fPrimers", "--forw",
		type=str,
		default=None,
		required=True,
		help="Path to forward primers fasta file")

	parser.add_argument("-rPrimers", "--rev",
		type=str,
		default=None,
		required=True,
		help="Path to reverse primers fasta file")

	parser.add_argument("-maxIndel", "--indel",
		type=int,
		default=None,
		required=True,
		help="Maximum indel (in nt) for read mapping")

	parser.add_argument("-minRatio", "--ratio",
		type=float,
		default=None,
		required=True,
		help="Minimum percent of read length to permit mapping")

	parser.add_argument("-o", "--out_dir",
		type=str,
		default=None,
		required=True,
		help="Directory to save output")


	args = parser.parse_args()

	numTargets = args.numTargs
	refSeqs = args.ref_seqs
	read1 = args.read1
	read2 = args.read2
	primF = args.forw
	primR = args.rev
	maxIndel = args.indel
	minRatio = args.ratio
	out_dir = args.out_dir

	# create the output directory if needed
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	checkVersion()
# 	makeOutDirs(out_dir)
# 	ref_fasta = getRef(ref_dir)
# 	sample = align(ref_fasta, pair1, pair2, out_dir)
# 	varCall(ref_fasta, sample, baseQual, out_dir)
# 
