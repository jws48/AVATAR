#!/usr/bin/env python

import os
import subprocess
import argparse
import glob
import sys

"""
This script is called as a single sample submission or Slurm array job, and takes in
the directory of the reference sequence,  forward read, reverse read, minimum base 
quality, and the output directory. If users do not submit a minumum base quality 
used for filtering variants, a phred score of 20 will be applied.

Sorted BAM and VCF files are output into the specified directory.

*************************************************************************************
Forward and reverse reads must be stored in separate directories.
*************************************************************************************

command: python variantArray.py -ref <ref_dir> -p1 <pair1> -p2 <pair2> -bq <baseQual> -o <out_dir>
"""

#Functions--------------------------------------------------------------------------
def checkVersion():
	if sys.version_info[0] < 3:
		raise Exception("Must be using Python 3")
		print("\nType: module load python\n")

def makeOutDirs(topOut):
	sam = os.path.join(topOut, "sam")
	bam = os.path.join(topOut, "bam")
	sortBam = os.path.join(topOut, "sort_bam")
	vcf = os.path.join(topOut, "vcf")
	vcf_unfilt = os.path.join(topOut, vcf, "unfiltered")
	vcf_filt = os.path.join(topOut, vcf, "filtered")
	if not os.path.exists(sam):
        	os.makedirs(sam)
	if not os.path.exists(bam):
		os.makedirs(bam)
	if not os.path.exists(sortBam):
		os.makedirs(sortBam)
	if not os.path.exists(vcf):
		os.makedirs(vcf)
	if not os.path.exists(vcf_unfilt):
		os.makedirs(vcf_unfilt)
	if not os.path.exists(vcf_filt):
		os.makedirs(vcf_filt)
		

def getRef(refDir):
	i = 0
	for f in os.listdir(refDir):
		if f.endswith((".fasta", ".fa")):
			fasta = os.path.join(refDir, f)
			i+=1
		elif f.endswith((".amb", ".ann", ".bwt", ".fai", ".pac", ".sa")):
			i+=1
		
	if i < 7:
		print("\n***Indexing reference sequence***\n")
		os.system("bwa index {}".format(fasta))
		os.system("samtools faidx {}".format(fasta))
	else:
		print("\n***Reference sequence already indexed***\n")
	return(fasta)


def align(refSeq, read1, read2, out):

	filename = os.path.basename(read1)
	sample = filename.split('_')[0]
	filename2 = os.path.basename(read2)
	sample2 = filename2.split('_')[0]

	assert sample == sample2
	
	os.system("bwa mem {} {} {} > {}/sam/{}.sam".format(refSeq, read1, read2, out, sample))
	os.system("samtools view -b {}/sam/{}.sam > {}/bam/{}.bam".format(out, sample, out, sample))
	os.system("rm {}/sam/{}.sam".format(out, sample))
	os.system("samtools sort {}/bam/{}.bam -o {}/sort_bam/{}_sort.bam".format(out,sample,out,sample))
	os.system("rm {}/bam/{}.bam".format(out, sample))

	return(sample)

def varCall(refSeq, sample, bq, out):
	anno = "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR"

	print("\nCalling variants\n")
	os.system("bcftools mpileup --annotate {} -Q {} -d 300000 -f {} {}/sort_bam/{}_sort.bam \
		| bcftools call -m --ploidy 1 --skip-variants indels > \
		{}/vcf/unfiltered/{}.vcf".format(anno, bq, refSeq, out, sample, out, sample))
	os.system("bcftools filter -g 1 {}/vcf/unfiltered/{}.vcf > {}/vcf/filtered/{}.vcf".format(out, sample, out, sample))
	print("Filtered VCF file in {}/vcf/filtered".format(out))

#Main-------------------------------------------------------------------------------
if __name__ == "__main__":

	description = "Calls variants on amplicon sequence data " \

	parser = argparse.ArgumentParser(prog="Variant caller", description=description)
	# user-defined parameters
	parser.add_argument("-ref", "--ref_dir",
		type=str,
		default=None,
		required=True,
		help="path to reference fasta file directory")

	parser.add_argument("-p1", "--pair1",
		type=str,
		default=None,
		required=True,
		help="path to forward read")

	parser.add_argument("-p2", "--pair2",
		type=str,
		default=None,
		required=True,
		help="path to reverse read")

	parser.add_argument("-bq", "--baseQual",
		type=int,
		default=20,
		required=False,
		help="the minimum base quality for calling a variant")

	parser.add_argument("-o", "--out_dir",
		type=str,
		default=None,
		required=True,
		help="directory to save output")


	args = parser.parse_args()

	ref_dir = args.ref_dir
	pair1 = args.pair1
	pair2 = args.pair2
	baseQual = args.baseQual
	out_dir = args.out_dir

	# create the output directory if needed
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	print("\nMinimum base quality set to {}\n".format(baseQual))
	checkVersion()
	makeOutDirs(out_dir)
	ref_fasta = getRef(ref_dir)
	sample = align(ref_fasta, pair1, pair2, out_dir)
	varCall(ref_fasta, sample, baseQual, out_dir)

