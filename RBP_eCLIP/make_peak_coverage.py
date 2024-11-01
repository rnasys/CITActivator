#! usr/bin/env python

"""
Author: Cikiy Wang
Date: 2020-8-31
E-mail: wangsiqi@picb.ac.cn
Description: extract single nucleotide coverage from IDR peak files
"""

"""
input file : 
    a. peak file :
		chr1    start  end  +/-
		... 	...   	...
"""

import argparse
import pandas as pd
import numpy as np
from sympy import *
import os
import sys
from multiprocessing import Process

def createHelp():
	"""
		Create the command line interface of the program.
	"""
	epilog_string = "Any question is welcome reported to wangsiqi@picb.ac.cn"
	description_string = 'The program is going to xtract single nucleotide coverage from IDR peak files'
	parser = argparse.ArgumentParser(description = description_string,epilog = epilog_string)
	parser.add_argument('-rbp', '--input RBP', dest = 'rbp', required=True, type = str,help = 'input RBP sample of eCLIP to analysis')
	parser.add_argument('-c', '--input celllines', dest = 'cell', required=True, type = str,help = 'input eCLIP data celllines')
	parser.add_argument('-r', '--input replicates', dest='rep', required=True, type=str,help='input replicates of eCLIP experiment to analysis (1 or 2)')
	parser.add_argument('-i_dir', '--input-dir', dest = 'indir', required=True, help = 'input peak bed files dir')
	parser.add_argument('-o_dir', '--output_dir', dest = 'outdir', required=True, help = 'output coverage files dir')
	opt = parser.parse_args()
	return opt

opt = createHelp()

RBP = opt.rbp
cellline = opt.cell
rep = opt.rep
peak_bed4_dir = opt.indir
peak_coverage_dir = opt.outdir

bed4_file = peak_bed4_dir + RBP + "_" + cellline + "_" + rep + ".bed4"
peak_coverage_file = peak_coverage_dir + RBP + "_" + cellline + "_" + rep + ".coverage"

bed = pd.read_csv(bed4_file, sep = "\t", names = ["chr","start","end","strand"])

#os.system("rm " + peak_coverage_file)
outfile = open(peak_coverage_file,'a')

outfile.write("chr\tpos\tcoverage\n")

for index,row in bed.iterrows():
    chr = row["chr"]
    start = row["start"]
    end = row["end"]
    for i in range(int(start),int(end)+1):
        outfile.write(str(chr) + "\t" + str(i) + "\t" + "1" + "\n")


