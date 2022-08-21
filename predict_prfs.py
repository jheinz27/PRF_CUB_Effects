import numpy as np
import argparse
import os
import math

parser = argparse.ArgumentParser(description='Check coverages in given region')
parser.add_argument('-w', required=True, help= 'Enter desired Window Length in NT, must be multiple of 3' )  
parser.add_argument('-fi', required=True, help= 'Fasta File of Sequences') 
parser.add_argument('-o', required=False, default = 'output.txt',help= 'output file') 
parser.add_argument('-b', required=False, default = '15', help= 'Enter desired buffer in NT, must be multiple of 3' )  
parser.add_argument('-p', required=False, default = '0.05', help= 'percentile' )  

args = parser.parse_args()
window = int(args.w) 
fastaFile = args.fi
output = args.o
buffer = int(args.b)
percentile = float(args.p)
#initialize type to int!!! 


#get aminoacid conversions from https://raw.githubusercontent.com/zhanxw/anno/master/codon.txt 
dict_aa = {}
with open('codon.txt', 'r') as f: 
    for line in f: 
        if not line.startswith('#'):
            s = line.strip().split('\t')
           
            dict_aa[s[0]] = s[2] 
            
            
#read in relative adaptive values from https://www.biologicscorp.com/tools/CAICalculator/ 
rel_weights = {} 
with open('yeast_rel_adapt_adj.txt', 'r') as f: 
    for line in f:  
        s = line.strip().split('\t')
        rel_weights[s[0]] = float(s[4])  
        
#define stop codons: 
stop_codons = []
for k in dict_aa.keys(): 
    if dict_aa[k] == 'O':
        stop_codons.append(k)

#read in elongation rates
elongate_dic = {}
with open('elongationrate.csv', 'r') as f: 
    for line in f:
        s = line.strip().split(';')
        elongate_dic[s[2]] = float(s[3])
       
#count codon frequencies
def get_codon_occurences(sequence): 
    codon_dict = {}

    
    #iterate through sequence in 3's to grab new codon
    i = 0
    while (i+3) < len(sequence):
        
        codon = sequence[i:i + 3]
       
        #increment frequency of codon appearances by one
      
        if codon not in codon_dict.keys():
            codon_dict[codon] = 1
        else:
            codon_dict[codon] += 1
        i += 3
        
    return codon_dict
    
#clean for early stop codons
def adj_for_stop(sequence): 
    i = 0
    while (i+3) < len(sequence):
        
        codon = sequence[i:i + 3]
       
        #increment frequency of codon appearances by one
        if codon in stop_codons: 
            return sequence[:i]
        
        i += 3 
    return sequence
    
    
#method to calc CAI
def calc_cai_cost(genomic_seq): 
    codon_freq = get_codon_occurences(genomic_seq)
    precai = 0 
    cost = 0 
    for k in codon_freq.keys():
        precai += codon_freq[k]*np.log(rel_weights[k])
        cost += codon_freq[k] * elongate_dic[k] 

    return np.exp(precai / (len(genomic_seq) /3)), cost

#method to combine potentials into a region with a buffer
def get_regions(pots):
    sort = sorted(pots)
    all_regions = []
    
    curset = []
    
    for i in range(1,len(sort)):
       
        cur = sort[i]
        prev = sort[i-1]
        curset.append(prev)
        if cur - prev <= buffer:
            
            curset.append(cur)
          
        else: 
            if len(curset) > 1 and (np.max(curset) - np.min(curset)) >= 9: 
                all_regions.append([int(np.min(curset)- (buffer/3)), int(np.max(curset) + (buffer/3))])
                #all_regions.append([np.min(curset), np.max(curset)])
            curset = []
    if len(curset) > 1 and (np.max(curset) - np.min(curset)) >= 9:
        all_regions.append([int(np.min(curset)- (buffer/3)), int(np.max(curset) + (buffer/3))])
        #all_regions.append([np.min(curset), np.max(curset)]) 
            
    return all_regions
            

        

#read in sequences
sequences = {} 
with open(fastaFile, 'r') as g: 
    for line in g: 
        l = line.strip()
        if l.startswith('>'):
            info = l.split(':')
            gene_name = info[-1]
        else: 
            raw_seq = l.upper()
            sequences[gene_name] = adj_for_stop(raw_seq)
           


#calculate cai scores for every possible position in every mRNA
potentials = {}
regions = {}
dcai = {} 
for k in sequences.keys(): 
    seq = sequences[k]
    pos_cai = {}
    cais = []
    pos = []
    #shifting windows
    for i in range(window, len(seq) - window, 3):
       
        before = seq[(i - window):i]
        after = seq[i +6 : (i + window +6) ]


        cai_before, cost_before = calc_cai_cost(before.upper())
        cai_after,cost_after  = calc_cai_cost(after.upper())
        delta_cai = cai_after - cai_before
        pos.append(i)
        cais.append(delta_cai)
    
    #percentile of interest
    top = math.ceil(percentile * len(cais))
    
    #to get top percentile
    maxes= np.argpartition(np.array(cais), -top)[-top:]
    #to get bottom percentile 
    #maxes= np.argpartition(np.array(cais), top)[:top]
    
    pos = np.array(pos)
    cais = np.array(cais)
    
    #get info about these maximum sites
    if len(maxes) > 0:
        regs = get_regions(pos[maxes])
        
    
        #potentials[k] = regs
        potentials[k] = sorted(pos[maxes])
        regions[k] = regs
        dcai[k] = sorted(cais[maxes])


    
#read in prf sites
id2prfloc = {}
with open('prfdb_S228C_all_nupack_unique_DNA_filtered_for_all.tsv','r') as f: 
    for line in f: 
        s = line.strip().split()
        info = int(s[1])
        if s[0] in id2prfloc.keys(): 
            temp = id2prfloc[s[0]]
            temp.append(info)
            id2prfloc[s[0]] = temp
        else: 
            id2prfloc[s[0]] = [info]

        
#write out regions of interest
with open(output, 'w') as out:
    out.write('#' + str(args) + '\n')
    for k in potentials.keys(): 
        if k in id2prfloc.keys(): 
            all_prfs = sorted(id2prfloc[k])
            all_regs = regions[k]
            all_cais = dcai[k]
            
            cur = potentials[k] 
            
            pots = str(cur[0])
            for i in range(1, len(cur)):
                pots += ',' + str(cur[i])
                
            prfs = str(all_prfs[0])
            for i in range(1, len(all_prfs)):
                prfs += ',' + str(all_prfs[i])
            
        
        
            
            if len(all_regs) == 0:
                regString = 'n/a'
            else:
                regString = str(all_regs[0][0]) + '-' + str(all_regs[0][1])
            for i in range(1,len(all_regs)):
                regString += ',' + str(all_regs[i][0]) + '-' + str(all_regs[i][1])
                
            caiString =  str(round(all_cais[0],4))
            for i in range(1, len(all_cais)):
                caiString += ',' + str(round(all_cais[i], 4))
                
            out.write(k + '\t' + pots + '\t' + caiString + '\t' + regString + '\t' + prfs + '\n') 
                
        
    
 
    
    
        

        
        
        
        
        
            
