#!/usr/bin/env python
# coding: utf-8

# #### You will need to read the sequences from the file and concatenate multiple lines into a single sequence for each identifier.
# For the first part of this problem, you will write code to compute codon usage statistics for the “surface glycoprotein" (spike protein) of the viral coding sequences in covid19codingClean.fasta. The data file is a text file containing sequences in fasta format: a header line starting with the \>" character that contains a description of the forthcoming sequence, and then the sequence itself, spread out over many lines. 

from Bio import SeqIO
s_protein_seq = []
seq_name = []

for seq_record in SeqIO.parse('covid19codingClean.fasta', "fasta"):
    genename = seq_record.description.split('|')[1].split('[')[0]
    if genename == "surface glycoprotein ":
        seq_name.append(seq_record.name)
        s_protein_seq.append(str(seq_record.seq))

codon_table = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
           "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
           "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
           "TTG" : "L", "CTG" : "L", "ATG" : "Start", "GTG" : "V",
           "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
           "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
           "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "TAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "TAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
           "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "TGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" }

CodonsDict = {  "TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0, "CTT": 0, 
        "CTC": 0, "CTA": 0, "CTG": 0, "ATT": 0, "ATC": 0, 
        "ATA": 0, "ATG": 0, "GTT": 0, "GTC": 0, "GTA": 0, 
       "GTG": 0, "TAT": 0, "TAC": 0, "TAA": 0, "TAG": 0, 
       "CAT": 0, "CAC": 0, "CAA": 0, "CAG": 0, "AAT": 0, 
       "AAC": 0, "AAA": 0, "AAG": 0, "GAT": 0, "GAC": 0, 
       "GAA": 0, "GAG": 0, "TCT": 0, "TCC": 0, "TCA": 0, 
       "TCG": 0, "CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0, 
       "ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0, "GCT": 0, 
       "GCC": 0, "GCA": 0, "GCG": 0, "TGT": 0, "TGC": 0, 
       "TGA": 0, "TGG": 0, "CGT": 0, "CGC": 0, "CGA": 0, 
       "CGG": 0, "AGT": 0, "AGC": 0, "AGA": 0, "AGG": 0, 
       "GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0} 

def count_codon(seq):
    dna =  seq
    self_dict = CodonsDict.copy()
    for i in range(0, len(dna)-(3+len(dna)%3), 3):
        codon = dna[i:i + 3]
        if codon in codon_table:
            self_dict[codon] +=1
            if codon_table[codon] == 'STOP':
                break
    
    return self_dict

import collections
col = ['accnum']
codon_dict = collections.OrderedDict(sorted(count_codon(s_protein_seq[0]).items()))
for k,v in codon_dict.items():
    col.append(k)

output = open('FPp3a-output','a')
tmp = ','.join(col)
output.write(tmp+'\n')
output.close()

output = open('FPp3a-output','a')
for i in range(0, len(s_protein_seq)):
    res=[]
    res.append(seq_name[i])
    codon_dict = collections.OrderedDict(sorted(count_codon(s_protein_seq[i]).items()))
    for k,v in codon_dict.items():
        res.append(v)
    output.write(str(res).split('[')[1].split(']')[0]+'\n')
output.close()


with open('FPp3a-output','r') as f:
    for num, line in enumerate(f):
        if num == 1:
            first_save = line.strip('\n')
            break
observe = first_save.split(',')[1:]
observe = [int(x) for x in observe]

n = sum(observe)
ref = {"A": 0.290, "C": 0.189, "G": 0.180, "T": 0.330}
exp = []
for i in range(1,len(col)):
    prob = 1.0
    tmp = col[i]
    for j in range(0,3):
        if tmp[j] in ref:
            prob *= ref.get(tmp[j])
    exp.append(prob*n)

import numpy as np
obs = np.array(observe)
chi_suqare = ['NC_045512.2']
for i in range(0,len(exp)):
    X = (np.square(obs[i] - exp[i]))/exp[i]
    chi_suqare.append(X)

col[0] = 'chi_suqare'
with open('FPp3b','a') as f:
    tmp = ','.join(col)
    f.write(tmp+'\n')
    f.write(str(chi_suqare).split('[')[1].split(']')[0]+'\n')

chi_sort = np.argsort(chi_suqare[1:])[::-1]
last_codon = []
last_idx = []
for i in range(0,len(chi_sort)):
    if obs[chi_sort[i]]<=exp[chi_sort[i]]:
        last_codon.append(col[chi_sort[i]+1])
        last_idx.append(chi_sort[i])

with open('FPp3b-rare','a') as f:
    for i in range(0,5):
        f.write(last_codon[i])
        f.write('\t'+str(round(chi_suqare[last_idx[i]]))+'\n')