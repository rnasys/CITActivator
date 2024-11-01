import pandas as pd
import sys
import re

#
# usage: python overlap_peaks.py sample region(150)
# output format: 
#

def intervalIntersection(chr,A,B):
    res = ""
    i = 0
    overlap = ""
    info_A = ""
    info_B = ""
    while i < int(A.iloc[:,0].size):
        j = 0
        while j < int(B.iloc[:, 0].size):
            #print(A["start"].iloc[1])
            low = max(int(A["start"].iloc[i]), int(B["start"].iloc[j]))
            high = min(int(A["end"].iloc[i]), int(B["end"].iloc[j]))
            #print(A["strand"].iloc[i])
            #print(B["strand"].iloc[j])
            if low <= high and (A["strand"].iloc[i] == B["strand"].iloc[j]):
                #print("low")
                info_A = str(chr) + ":" + str(A["start"].iloc[i]) + "-" + str(A["end"].iloc[i])
                info_B = str(chr) + ":" + str(B["start"].iloc[j]) + "-" + str(B["end"].iloc[j])
                overlap = str(chr) + ":" + str(low) + "-" + str(high)
                overlap_len = str(high-low)
                res = res + info_A + "\t" + info_B + "\t" + overlap + "\t" + overlap_len + "\n" 
            j += 1
        i += 1
    return res

def peak2bed_1(peak): # format: chr start end strand
    re_peak = peak[["chr","start","end","strand"]]
    re_peak_chr = re_peak["chr"].unique()
    return re_peak,re_peak_chr

def peak2bed_2(peak): # format: chr start end strand
    re_peak_pos = peak["Location (hg19)"].str.split(";",2,True)
    re_peak_new = pd.DataFrame()
    re_peak_new["pos"] = re_peak_pos.iloc[:,0].tolist()+re_peak_pos.iloc[:,1].dropna().tolist()
    re_peak_strand = peak["Strand (hg19)"].str.split(";",2,True)
    re_peak_new["strand"] = re_peak_strand.iloc[:,0].tolist()+re_peak_strand.iloc[:,1].dropna().tolist()
    
    re_peak_new = re_peak_new.dropna()
    re_peak = re_peak_new["pos"].str.split("[:-]",3,True)
    re_peak = re_peak.dropna()
    
    #re_peak = re_peak.drop([3],axis=1)
    re_peak.columns = ["chr","start","end"]
    re_peak["start"] = [int(c.replace(',',''))-int(region) for c in re_peak["start"]]
    re_peak["end"] = [int(c.replace(',',''))+int(region) for c in re_peak["end"]]
    re_peak["strand"] = re_peak_new["strand"]
    re_peak_chr = re_peak["chr"].unique()
    return re_peak,re_peak_chr

# you need to change your work directory
sample = sys.argv[1]
region = sys.argv[2]
peak_file_1 = ""+sample+"_HepG2_rep12.bed6"
peak_file_2 = "Human_IRES_Info.txt"
work_dir = ""

# main
peak_1 = pd.read_csv(peak_file_1,sep = "\t",encoding="unicode_escape",names=["chr","start","end","pos","ids","strand"],index_col=False)
peak_2 = pd.read_csv(peak_file_2,sep = "\t",encoding="unicode_escape")

peak_1_bed = peak2bed_1(peak_1)[0]
peak_2_bed = peak2bed_2(peak_2)[0]
a = peak2bed_1(peak_1)[1]
b = peak2bed_2(peak_2)[1]
chr_overlap = set(a).intersection(set(b))

#print(chr_overlap)

results = "peak\tires\toverlap\toverlap_length\n"
for chr in chr_overlap:
    peak_1_bed_chr = peak_1_bed[peak_1_bed["chr"] == chr]
    peak_2_bed_chr = peak_2_bed[peak_2_bed["chr"] == chr]
    #print(sig_peak_DiffBind.iloc[:, 0].size)
    #print(sig_peak_tcga.iloc[:,0].size)
    print(chr + " START")
    re = intervalIntersection(chr, peak_1_bed_chr, peak_2_bed_chr)
    results = results + re
    print(chr + " DONE")

file = open(work_dir + sample +"_peak_ires_overlap.txt", "w")
file.write(results)
file.close()





