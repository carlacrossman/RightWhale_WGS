#!/usr/bin/env python3

from __future__ import print_function

import inspect
import argparse
import sys
import datetime
import os
import gzip
import vcfpy
import random


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def parseKeep(keep):
	keepList = []
	file=open(keep, "r+")
	for line in file.readlines():
		line = line.rstrip()
		keepList.append(line)
	return keepList

def readVCF(vcfFile, minFreq, keep, randomAllele, missdata):

	# Open file, this will read in the header
	reader = vcfpy.Reader.from_path(vcfFile)
	reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'AA'), ('Number', '1'), ('Type', 'String'), ('Description', 'Ancestral Allele')]))
	# writer = vcfpy.Writer.from_path('/dev/stdout', reader.header)

	keepListIndex = []
	if keep:
		keepList = parseKeep(keep)
		for ind in keepList:
			[keepListIndex.append(i) for i,x in enumerate(reader.header.samples.names) if x == ind]

	for record in reader:
		if not record.is_snv():
			continue

		# get the allele present in keeplistindex
		variant = {}
		printTest = "FALSE"
		bases = []
		none = 0
		k = 0 
		variant["aa_miss"] = 0
		plodity = 0

		for call in record.calls:
			if k in keepListIndex:
				if call.gt_alleles is None: # for haploid vcf none is reported instead of 
					variant["aa_miss"] += 1
					plodity = 1
				else :
					if len(call.gt_alleles) == 1: # for haploid vcf
						plodity = 1
						if call.gt_alleles[0] is None: # this case should not happen, here just in case
							variant["aa_miss"] += 1
						else:
							bases.append(call.gt_bases[0])

					elif len(call.gt_alleles) == 2: # for diploid vcf
						plodity = 2
						if call.gt_alleles[0] is None: # this case should not happen, here just in case
							variant["aa_miss"] += 1
						else:
							bases.append(call.gt_bases[0])

						if call.gt_alleles[1] is None: # this case should not happen, here just in case
							variant["aa_miss"] += 1
						else:
							bases.append(call.gt_bases[1])
			k += 1

		variant["aa_miss_proportion"] = variant["aa_miss"]/(variant["aa_miss"]+len(bases))
		nb = []
		nt = []
		for b in set(bases):
			nb.append(bases.count(b))
			nt.append(b)

		index =[]
		variant["aa_freq"] = 0
		variant["random_choice"] = "TRUE"
		[index.append(i) for i,x in enumerate(nb) if x == max(nb)]
		if len(index) == 1: # one single major variant
			if nb[index[0]]/len(bases)>minFreq:
				printTest = "TRUE"
				aa = nt[index[0]]
			variant["aa"] = nt[index[0]]
			variant["aa_freq"] = nb[index[0]]/len(bases)
			variant["random_choice"] = "FALSE"
		elif len(index) > 1: # equality among major variants
			r = random.choice(index)
			variant["aa"] = nt[r]
			variant["aa_freq"] = nb[r]/len(bases)
			variant["random_choice"] = "TRUE"
			if nb[r]/len(bases)>minFreq:
				if randomAllele:
					printTest = "TRUE"
					aa = nt[r]
				else:
					printTest = "FALSE"

		if record.REF in ["A", "T", "C", "G"] \
		and variant["aa_miss_proportion"] <= (1-missdata) \
		and variant["aa_freq"] >= minFreq \
		and variant["random_choice"] == "FALSE" :
			if variant["aa"] != record.REF:
				reorganizeVCF(record, variant["aa"], plodity)
			else:
				writeManuallyVCF(record, plodity)


def writeManuallyVCF(record, plodity):
	altList = []
	altList += [alt.value for alt in record.ALT]

	# for vcf having SY tag in info collum 
	info = "."
	if "SY" in record.INFO:
		info = str("SY="+record.INFO["SY"])
	if "4" in record.INFO:
		info += ";"+str("4="+record.INFO["4"])

	if plodity == 2:
		printList = [record.CHROM, str(record.POS), ".", record.REF, str.join(',', altList), str(record.QUAL), ".", str(info), "GT"]
		printList += [call.data.get('GT') for call in record.calls]
		print(str.join('\t', printList))

	elif plodity == 1:
		printList = [record.CHROM, str(record.POS), ".", record.REF, str.join(',', altList), str(record.QUAL), ".", str(info), "GT"]
		for call in record.calls:
			if call.gt_alleles is None:
				printList.append(".")
			else:
				printList.append(call.data.get('GT'))
		print(str.join('\t', printList))

def reorganizeVCF(record, aa, plodity):
	var = []
	var.append(aa)
	var.append(record.REF)

	for alt in record.ALT:
		if alt.value != aa:
			var.append(alt.value)
	alt = []
	for b in range(1,len(var)):
		alt.append(var[b])

	# for vcf having SY tag in info collum 
	info = "."
	if "SY" in record.INFO:
		info = str("SY="+record.INFO["SY"])
	if "4" in record.INFO:
		info += ";"+str("4="+record.INFO["4"])

	printList = [record.CHROM, str(record.POS), ".", aa, str.join(',', alt), str(record.QUAL), ".", str(info), "GT"]
	gtList = []

	if plodity == 2:
		for call in record.calls:
			# print(call.gt_bases[0])
			if call.gt_bases[0] is not None:
				gt = str(var.index(call.gt_bases[0]))+"/"+str(var.index(call.gt_bases[1]))
				gtList.append(gt)
				# print(gt)
			else:
				gt = "./."
				gtList.append(gt)
				# print(gt)
		printList += gtList
		print(str.join('\t', printList))

	elif plodity == 1:
		for call in record.calls:
		# print(call.gt_bases[0])
			if call.gt_alleles is None:
				gt = "."
				gtList.append(gt)

			else:
				gt = str(var.index(call.gt_bases[0]))
				gtList.append(gt)
				# print(gt)
		printList += gtList
		print(str.join('\t', printList))

if __name__=="__main__": # main
	pgrmName = os.path.basename(__file__)
	parser=argparse.ArgumentParser(
		prog=pgrmName,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='''
Progam is is made to polarize a vcf file.
To polarize a vcf file the program will take the most frequent allele
present across all individuals or present across a selected list of individuals.''',
		epilog="""...""")

	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	parser.add_argument("--vcf", dest="vcf", metavar="FILE", required=True,
		help="Defines the VCF file to be processed (either .vcf or .vcf.gz).")

	parser.add_argument("--keep", dest="keep", metavar="FILE", required=False,
		help="Provide files containing a list of individuals in which the ancestral allele will be determined.")

	parser.add_argument("--miss", dest="missdata", type=float, default=1,
		help="Proportion of AA missing data allow, 1=no missing data allowed [0:1] default(1)")

	parser.add_argument("--minFreq", dest="minFreq", type=float, default=0,
		help="Allele minimun frequency to be considered as AA [0:1] default(0)")

	parser.add_argument("-r", dest="random", action='store_true',
		help="If equal proportion shared by multiple allele, choose randomly one allele over the other.")

	## if no args then return help message
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)

	args = parser.parse_args()
	# readVCF(args.vcf, args.missdata, args.minFreq, args.keep, args.random)
	readVCF(args.vcf, args.minFreq, args.keep, args.random, args.missdata)