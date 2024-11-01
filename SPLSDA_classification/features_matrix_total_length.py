import pandas as pd
import numpy as np

'''
    usage : python features_matrix_total_length.py
    input single-peptide/di-peptides/tri-peptides counts in RBP_domain_aa_IDR 
    output matrix : RBP(110) x Features(4963)
'''

dir_kmer_disorder_stat = ""
dir_single_disorder_stat = ""
dir_data = ""
kmer2 = pd.read_csv(dir_kmer_disorder_stat + "RBP_domain_aa_IDR_total_kmer2.stat",sep = "\t")
kmer3 = pd.read_csv(dir_kmer_disorder_stat + "RBP_domain_aa_IDR_total_kmer3.stat",sep = "\t")
kmer2_bg = pd.read_csv(dir_kmer_disorder_stat + "1542_RBP_IDR_total_kmer2.stat",sep = "\t")
kmer3_bg = pd.read_csv(dir_kmer_disorder_stat + "1542_RBP_IDR_total_kmer3.stat",sep = "\t")
sg_peptide = pd.read_csv(dir_single_disorder_stat + "RBP_domain_total.stat",sep = "\t")
sg_peptide_bg = pd.read_csv(dir_single_disorder_stat + "RBP1542_total.stat",sep = "\t")

sg_peptide.index = sg_peptide["Unnamed: 0"].values.tolist()
sg_peptide = sg_peptide.drop("Unnamed: 0",axis=1)
x = sg_peptide.div(sg_peptide.sum(axis=0), axis=1) # sum(axis=0) colsum

sg_peptide_bg.index = sg_peptide_bg["Unnamed: 0"].values.tolist()
sg_peptide_bg = sg_peptide_bg.drop("Unnamed: 0",axis=1)
temp = sg_peptide_bg.sum(axis=1)
bg_freq = temp/temp.sum(axis=0)
x_bg = x.div(bg_freq, axis=0)

kmer2.index = kmer2["Unnamed: 0"].values.tolist()
kmer2 = kmer2.drop("Unnamed: 0",axis=1)
y = kmer2.div(kmer2.sum(axis=0), axis=1)

kmer2_bg.index = kmer2_bg["Unnamed: 0"].values.tolist()
kmer2_bg = kmer2_bg.drop("Unnamed: 0",axis=1)
temp = kmer2_bg.sum(axis=1)
bg_freq = temp/temp.sum(axis=0)
y_bg = y.div(bg_freq, axis=0)
y_bg = y_bg.dropna()

kmer3.index = kmer3["Unnamed: 0"].values.tolist()
kmer3 = kmer3.drop("Unnamed: 0",axis=1)
z = kmer3.div(kmer3.sum(axis=0), axis=1)

kmer3_bg.index = kmer3_bg["Unnamed: 0"].values.tolist()
kmer3_bg = kmer3_bg.drop("Unnamed: 0",axis=1)
temp = kmer3_bg.sum(axis=1)
bg_freq = temp/temp.sum(axis=0)
z_bg = z.div(bg_freq, axis=0)
z_bg = z_bg.dropna()

freq_matrix = np.transpose(x.append(y).append(z))
#freq_matrix = np.transpose(x.append(y))
#freq_matrix = np.transpose(z)
#freq_matrix = np.transpose(x_bg.append(y_bg).append(z_bg))
freq_matrix = freq_matrix.dropna(axis=0)
#print(x)
#print(y)
#print(z)
print(freq_matrix)

#count_matrix.columns = count_matrix.iloc[[0]].values.tolist()
#count_matrix = count_matrix.drop(index=["Unnamed: 0"])
#freq_matrix = freq_matrix.drop(["DTF","EIQ"],axis = 1)
freq_matrix.to_csv(dir_data + "feature_matrix_RBPdomain_total_length.txt",sep = "\t")

