#! usr/bin/env python

"""
Author: Siqi Wang
Date: 2020-2-29
E-mail: wangsiqi@picb.ac.cn
Description: normalized depth on full transcripts based on eCLIP data
"""

"""
input file : 
	a. annotation file with comprehensive format : 
		name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds
	b. coverage file :
		chr1    142680  1
		chr1    238754  1
		chr1    334173  1
		... 	...   	...
	c. expression results file from rsem : 
		
	d. protein coding transcripts annotation file :

version 1.1:
	* step_size win_size
	get_high_exp_transcripts
	max
	win: 13:100:49
	norm = depth x trans_depth / (exon_len x total_trans / 1000)
	SMInput
	transcript_aver_depth

use command:
	python v3.1_single_nucleotide_depth_norm_single_threading.py \
		-depth_dir /dir/ \
		-rbp HNRNPUL1 \
		-rep 1 \
		-i_expr rsem.isoforms.results \
		-o /dir/ > HNRNPUL1_1.v3.1.log

"""
import argparse
import pandas as pd
import numpy as np
import time
from sympy import *
from multiprocessing import Process

def createHelp():
	"""
		Create the command line interface of the program.
	"""
	epilog_string = "Any question is welcome reported to wangsiqi@picb.ac.cn"
	description_string = 'The program is going to get normalized depth on full transcripts based on eCLIP data'
	parser = argparse.ArgumentParser(description = description_string,epilog = epilog_string)
	parser.add_argument('-rbp', '--input RBP', dest = 'rbp', required=True, type = str,help = 'input RBP sample of eCLIP to analysis')
	parser.add_argument('-rep', '--input replicates', dest = 'rep', required=True, type = str,help = 'input replicates of eCLIP experiment to analysis (1 or 2)')
	parser.add_argument('-i_ref', '--input-ref', dest = 'anno', default = '/picb/rnasys2/wangsiqi/Annotation/human/Gencode_v19/Gencode_v19_comprehensive.txt', help = 'input annotation file')
	parser.add_argument('-i_pc', '--input-pc', dest = 'pc', default = '/picb/rnasys2/wangsiqi/Annotation/human/Gencode_v19/gencode.v19.protein_coding_gene_transcript.bed', help = 'input annotation file')
	parser.add_argument('-i_expr', '--input-expr', dest = 'expr', required=True, help = 'input transcripts expression file from rsem')
	parser.add_argument('-depth_dir', '--depth dir', dest = 'depth_dir', required=True, help = 'input project files directory')
	#parser.add_argument('-cover_dir', '--coverage dir', dest = 'cover_dir', required=True, help = 'input coverage files directory')
	parser.add_argument('-w_UTR5', '--win-UTR5', dest = 'w_utr5', default = 13, type = int,help = 'input bin number in UTR5')
	parser.add_argument('-w_CDS', '--win-CDS', dest = 'w_cds', default = 100, type = int,help = 'input bin number in CDS')
	parser.add_argument('-w_UTR3', '--win-UTR3', dest = 'w_utr3', default = 49, type = int,help = 'input bin number in UTR3')
	#parser.add_argument('-s', '--strand', dest = 's', default = '+', type = str, help = 'strand (+ or -)')
	#parser.add_argument('-p', '--thread', dest = 'p', default = 10, type = int, help = 'thread number used in multiprocessing')
	parser.add_argument('-o', '--output dir', dest = 'Out_dir', required=True, type = str,help = 'output directory')
	opt = parser.parse_args()
	return opt

def get_max_transcripts(rsem_pc):
	"""
		Get max expressed protein-coding transcripts of each gene for calculation from rsem results 
	"""
	isoforms_0 = rsem_pc[rsem_pc["IsoPct"] > 0 ] #filter isoform exist expression
	isoforms_max = isoforms_0.groupby('gene_id').apply(lambda t: t[t.IsoPct==t.IsoPct.max()])[["transcript_id"]]
	pd_isoforms_max = pd.DataFrame(isoforms_max)
	print("\tget_max_transcripts done "+ time.asctime( time.localtime(time.time()) ))
	return pd_isoforms_max

def get_longest_transcripts(rsem_pc):
	"""
		Get longest protein-coding transcripts of each gene for calculation from rsem results 
	"""
	#isoforms_0 = rsem_pc[rsem_pc["IsoPct"] > 0 ] #filter isoform exist expression
	isoforms_0 = rsem_pc
	isoforms_max = isoforms_0.groupby('gene_id').apply(lambda t: t[t.effective_length==t.effective_length.max()])[["transcript_id"]]
	pd_isoforms_max = pd.DataFrame(isoforms_max)
	print("\tget_longest_transcripts done "+ time.asctime( time.localtime(time.time()) ))
	return pd_isoforms_max

def get_high_exp_transcripts(rsem_pc):
	"""
		Get highly expression protein-coding transcripts of each gene for calculation from rsem results （TPM>=1）
	"""
	isoforms_0 = rsem_pc[rsem_pc["TPM"] >= 1 ] #filter isoform exist expression
	isoforms_max = isoforms_0.groupby('gene_id').apply(lambda t: t[t.IsoPct==t.IsoPct.max()])[["transcript_id"]]
	pd_isoforms_max = pd.DataFrame(isoforms_max)
	print("\tget_max_transcripts done "+ time.asctime( time.localtime(time.time()) ))
	return pd_isoforms_max

def get_coverage_transcripts(df):
	"""
		Get protein-coding transcripts of each gene
	"""
	df_0 = df[df["coverage"] > 0 ] #filter 
	df_max = df_0.groupby('gene_id').apply(lambda t: t[t.coverage==t.coverage.max()])[["transcript_id"]]
	pd_df_max = pd.DataFrame(df_max)
	print("\tget_longest_transcripts done "+ time.asctime( time.localtime(time.time()) ))
	return pd_df_max

def append_exons(exon_s,exon_e,exon_count):
	"""
		append all exons positions from annotation file (*_comprehensive.txt)
		exons format:
			33546755,33547778,33549554,33557650,···,
	"""
	exon_pos = []
	s = exon_s.split(",")
	e = exon_e.split(",")
	for i in range(0,int(exon_count)):
		new = np.array(range(int(s[i]),int(e[i]) + 1))
		exon_pos = np.append(exon_pos,new)
	print("\tappend_exons done "+time.asctime( time.localtime(time.time()) ))
	return exon_pos

def region_pos(start,end,pos,strand):
	"""
		sperate positions in 5'UTR, 3'UTR, CDS regions
		x : near start
		y : near end
	"""
	if strand == "+" :
		#also including transcripts with no UTR, short UTR...
		x = pos[pos >= start].min()
		y = pos[pos >= end].min()
		xloc = np.where(pos == x)[0][0]
		yloc = np.where(pos == y)[0][0]
		UTR5_pos = pos[0 : xloc]
		CDS_pos = pos[xloc : yloc]
		UTR3_pos = pos[yloc : len(pos)]
	elif strand == "-":
		pos_minus = pos[::-1]
		x = pos_minus[pos_minus <= end].max()
		y = pos_minus[pos_minus <= start].max()
		xloc = np.where(pos_minus == x)[0][0]
		yloc = np.where(pos_minus == y)[0][0]
		UTR5_pos = pos_minus[0 : xloc]
		CDS_pos = pos_minus[xloc : yloc]
		UTR3_pos = pos_minus[yloc : len(pos_minus)]
	print("	PASS: " + transcript + " has been calculated " + strand + " : ")
	return UTR5_pos,CDS_pos,UTR3_pos

def get_bin(length,bin_count):
	"""
		Get bin numbers with different size
	"""
	m = int(length/bin_count)
	n = m + 1
	x = Symbol('x')
	y = Symbol('y')
	sol = solve([m*x + n*y - length , x + y - bin_count])
	return m,n,sol[x],sol[y]

def bin_aver_depths(pos,bin_count):
	"""
		Calculate average depth in each bin
	"""
	if len(pos) >= bin_count:
		aver_depths = [0]*bin_count
		cal_bin = get_bin(len(pos),bin_count)
		small_bin_size = cal_bin[0]
		large_bin_size = cal_bin[1]
		small_bin_num = cal_bin[2]
		large_bin_num = cal_bin[3]
		for i in range(0,small_bin_num):
			pos_in_bin = pos[i*small_bin_size : (i+1)*small_bin_size]
			depth_sum_bin = 0
			for x in pos_in_bin:
				if x in pos_cov:
					depth = coverage_chr_tx_pos_index.at[x,"coverage"]
				else:
					depth = 0
				depth_sum_bin = depth_sum_bin + depth
			aver_depths[i] = depth_sum_bin/small_bin_size
			#np.intersect1d[pos_in_bin,pos_cov]
		for j in range(0,large_bin_num):
			pos_in_bin = pos[ small_bin_num*small_bin_size + j*large_bin_size: small_bin_num*small_bin_size + (j+1)*large_bin_size]
			depth_sum_bin = 0
			for x in pos_in_bin:
				if x in pos_cov:
					depth = coverage_chr_tx_pos_index.at[x,"coverage"]
				else:
					depth = 0
				depth_sum_bin = depth_sum_bin + depth
			aver_depths[small_bin_num+j] = depth_sum_bin/large_bin_size
	else:
		aver_depths = [0]*bin_count
	print("\tbin_aver_depths done "+time.asctime( time.localtime(time.time()) ))
	return np.array(aver_depths)

def win_aver_depths(pos,win_count):
	"""
		Calculate average depth in each window
	"""
	region_len = len(pos)
	if region_len < win_count: # UTR5 length is too short to split
		aver_depths = [0]*win_count
		input_aver_depths = [0]*win_count
	#elif region_len < win_count*10: # No overlap among each window
	#	win_size = int(region_len/win_count)
	#	step_size = int(region_len/win_count)
	#	aver_depths = [0]*win_count
	else:
		step_size = int(region_len/(win_count+1))
		win_size = int(region_len-(step_size*(win_count-1)))
		aver_depths = [0]*win_count
		input_aver_depths = [0]*win_count
		for i in range(win_count):
			win_start_pos=i*step_size
			pos_in_win = pos[win_start_pos : win_start_pos+win_size]
			depth_sum_win = 0
			input_depth_sum_win = 0
			for x in pos_in_win:
				if x in pos_cov:
					depth = coverage_chr_tx_pos_index.at[x,"coverage"] 
				else:
					depth = 0
				if x in input_pos_cov:
					input_depth = input_coverage_chr_tx_pos_index.at[x,"coverage"].sum()
				else:
					input_depth = 0
				depth_sum_win = depth_sum_win + depth
				input_depth_sum_win = input_depth_sum_win + input_depth
			aver_depths[i] = depth_sum_win/win_size
			input_aver_depths[i] = input_depth_sum_win/win_size
	print("\twin_aver_depths done "+time.asctime( time.localtime(time.time()) ))
	return np.array(aver_depths),np.array(input_aver_depths)

def transcript_aver_depth(pos):
	trans_exon_mapping_count = 0
	for x in pos:
		if x in pos_cov:
			depth = coverage_chr_tx_pos_index.at[x,"coverage"]
		else:
			depth = 0
		trans_exon_mapping_count = trans_exon_mapping_count + depth
	print("\ttranscript_depth done "+time.asctime( time.localtime(time.time()) ))
	return trans_exon_mapping_count

def depth_norm(aver_depth,count):
	#norm_aver_depth = (aver_depth*trans_exon_mapping_count)/(exon_length/1000)
	norm_aver_depth = (aver_depth*count)/(exon_length/1000)
	print("\tdepth_norm done "+time.asctime( time.localtime(time.time()) ))
	return norm_aver_depth

#__main__
print("START at "+time.asctime( time.localtime(time.time()) ))

opt = createHelp()
#output_dir = ""
#ann_file = "Gencode_v19_comprehensive.txt"
#pc_gene_trans_file = "gencode.v19.protein_coding_gene_transcript.bed"
#isoform_exp_file = "rsem.isoforms.results"

RBP = opt.rbp
rep = opt.rep
sample = RBP + "_" + rep
SMInput = RBP + "_input_1"

ann = pd.read_csv(opt.anno,sep = "\t",usecols = [1,2,3,4,5,6,7,8,9,10])
ann_drop = ann.drop_duplicates(["name"], "first")
ann_index = ann_drop.set_index('name')
anno_trans = ann["name"].tolist()
isoform_exp = pd.read_csv(opt.expr,sep = "\t")
pc_gene_trans = pd.read_csv(opt.pc,sep = "\t",names = ["gene_id","transcript_id"])
site_coverage = pd.read_csv(opt.depth_dir + sample + ".plus_minus.coverage.bed",sep = "\t",names = ["chr","pos","coverage"])
input_coverage = pd.read_csv(opt.depth_dir + SMInput + ".plus_minus.coverage.bed",sep = "\t",names = ["chr","pos","coverage"])

"""
	split into 24 chromosomes : ex: df_plus_chr1
"""
classinformation = site_coverage["chr"].unique()
for temp_classinformation in classinformation:
	temp_data = site_coverage[site_coverage["chr"].isin([temp_classinformation])]
	exec("df_%s = temp_data"%temp_classinformation)

in_classinformation = input_coverage["chr"].unique()
for in_temp_classinformation in in_classinformation:
	temp_data = input_coverage[input_coverage["chr"].isin([in_temp_classinformation])]
	exec("input_%s = temp_data"%in_temp_classinformation)

total_mapping_counts = site_coverage["coverage"].sum()
input_total_mapping_counts = input_coverage["coverage"].sum()

"""
	calculate the longest protein coding isofoms and select transcripts for following depth calculation
"""
isoform_exp_mg = pd.merge(pc_gene_trans, isoform_exp, on=['transcript_id','gene_id'])
isoform_exp_max = get_high_exp_transcripts(isoform_exp_mg)
all_transcripts = isoform_exp_max["transcript_id"]
total = len(all_transcripts)

#array with 0 values
norm_depth_tx_sum = np.array([0 for _ in range(opt.w_utr5 + opt.w_cds + opt.w_utr3)])

n = 1
#total_exon_mapping_count = 0
#total_mapping_count = site_coverage["coverage"].sum()
#input_total_mapping_count = input_coverage["coverage"].sum()

for transcript in all_transcripts:
	#print(transcript + " START at "+time.asctime( time.localtime(time.time()) ) + " (" + str(n) + "/" + str(total) +")")
	#print(anno_trans)
	if transcript in anno_trans:
		print(transcript + " START at "+time.asctime( time.localtime(time.time()) ) + " (" + str(n) + "/" + str(total) +")")
	else:
		continue
	"""
		information in annotation file
	"""
	chromsome = ann_index.at[transcript,"chrom"]
	txstart = int(ann_index.at[transcript,"txStart"])
	txend = int(ann_index.at[transcript,"txEnd"])
	tx_length = txend - txstart + 1
	strand = ann_index.at[transcript,"strand"]
	cdsstart = int(ann_index.at[transcript,"cdsStart"])
	cdsend = int(ann_index.at[transcript,"cdsEnd"])
	coverage_chr = eval("df_" + chromsome)
	coverage_chr_tx = coverage_chr[(coverage_chr["pos"] <= txend) & (coverage_chr["pos"] >= txstart)]
	coverage_chr_tx_pos_index = coverage_chr_tx.set_index('pos')
	input_coverge_chr = eval("input_" + chromsome)
	input_coverage_chr_tx = input_coverge_chr[(input_coverge_chr["pos"] <= txend) & (input_coverge_chr["pos"] >= txstart)]
	input_coverage_chr_tx_pos_index = input_coverage_chr_tx.set_index('pos')
	pos_cov = coverage_chr_tx.iloc[:,1].values # type = array
	input_pos_cov = input_coverage_chr_tx.iloc[:,1].values # type = array
	trans_mapping_counts = coverage_chr_tx["coverage"].sum()
	input_trans_mapping_counts = input_coverage_chr_tx["coverage"].sum()
	"""
		append exon positions --- sorted as increase
	"""
	#all_exons_positions = np.sort(append_exons(ann_index.at[transcript,"exonStarts"],ann_index.at[transcript,"exonEnds"],ann_index.at[transcript,"exonCount"])).astype(np.int64)
	all_exons_positions = append_exons(ann_index.at[transcript,"exonStarts"],ann_index.at[transcript,"exonEnds"],ann_index.at[transcript,"exonCount"]).astype(np.int64)
	#trans_exon_mapping_count = transcript_depth(all_exons_positions)
	#total_exon_mapping_count = total_exon_mapping_count + trans_exon_mapping_count
	exon_length = len(all_exons_positions)
	regions = region_pos(cdsstart,cdsend,all_exons_positions,strand)
	UTR_5 = regions[0]
	CDS = regions[1]
	UTR_3 = regions[2]
	#if cdsstart != cdsend:
	#else:
	#	print("CDS: " + transcript + "has no CDS regions")
	#	continue
	# calculate depth around start and stop codon in 3 regions
	UTR5_aver_depths = win_aver_depths(UTR_5,opt.w_utr5)
	CDS_aver_depths = win_aver_depths(CDS,opt.w_cds)
	UTR3_aver_depths = win_aver_depths(UTR_3,opt.w_utr3)
	# normalize depth around start and stop codon in 3 regions
	#norm_UTR5_aver_depths = depth_norm(UTR5_aver_depths)
	#norm_CDS_aver_depths = depth_norm(CDS_aver_depths)
	#norm_UTR3_aver_depths = depth_norm(UTR3_aver_depths)
	#append regions
	#norm_depth_tx = np.concatenate((norm_UTR5_aver_depths,norm_CDS_aver_depths,norm_UTR3_aver_depths),axis=0)
	depth_tx = np.concatenate((UTR5_aver_depths[0],CDS_aver_depths[0],UTR3_aver_depths[0]),axis=0)
	input_depth_tx = np.concatenate((UTR5_aver_depths[1],CDS_aver_depths[1],UTR3_aver_depths[1]),axis=0)
	norm_depth_tx = depth_norm(depth_tx,trans_mapping_counts)
	input_norm_depth_tx = depth_norm(depth_tx,input_trans_mapping_counts)
	relative_norm_depth = (norm_depth_tx/total_mapping_counts) - (input_norm_depth_tx/input_total_mapping_counts)
	relative_norm_depth[relative_norm_depth < 0] = 0
	#sum of all transcripts' normalized depth
	norm_depth_tx_sum = norm_depth_tx_sum + relative_norm_depth
	n = n + 1

print("DONE: All transcripts done!")

norm_depth_tx_sum_pd = pd.DataFrame(norm_depth_tx_sum)
norm_depth_tx_sum_pd.to_csv(path_or_buf = opt.Out_dir + sample + ".sn.norm.depth",sep = "\t")

print("DONE at "+time.asctime( time.localtime(time.time()) ))
