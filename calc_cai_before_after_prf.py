import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description='Check coverages in given region')
parser.add_argument('-w', required=True, help= 'Enter desired Window Length in NT, must be multiple of 3' )  
args = parser.parse_args()
window = int(args.w) 




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
        elongate_dic[s[2]] = 1/(float(s[3]) / 10)
       
        
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



sequences = {} 
#with open('yeast_gencode_all_exons_pcs_IDs_1L_filt.fa', 'r') as g: 
with open('expression_samsies.fa', 'r') as g: 
    for line in g: 
        if not line.startswith('-'):
            l = line.strip()
            if l.startswith('>'):
                info = l.split(':')
                gene_name = info[-1]
            else: 
                raw_seq = l.upper()
                sequences[gene_name] = adj_for_stop(raw_seq)
            
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
            
#window = 150



outfile = 'all_prfdb_cai_window'+ str(window)+'.tsv'

gap = 6
with open(outfile, 'w') as out:
    for k in id2prfloc.keys(): 
        if k not in sequences.keys():
            continue
            
        prfs = id2prfloc[k]
        seq = sequences[k]
        
        
        for p in prfs:
            
            curprf = p 
            if curprf <= window + gap or curprf >= (len(seq) - window - gap): 
                continue
            to_use = curprf - (curprf % 3) 
            
            before = seq[(curprf - window - gap ):curprf - gap ]
            after = seq[(curprf + gap): (curprf + window +gap)]
        

            cai_before, cost_before = calc_cai_cost(before.upper())
            cai_after,cost_after  = calc_cai_cost(after.upper())
            out.write(k + '\t' + str(p) + '\t' + str(len(seq)) +  '\t' + str(cai_before) + '\t' + str(cai_after) + '\t'+ str(cai_after-cai_before) +  '\t' + str(cost_before) + '\t' +str(cost_after) + '\t'+ str(round((cost_after - cost_before),4))+ '\n' )   
  


all_gene_ids = list(sequences.keys())

           
outPrf = 'all_random_locs_window' + str(window)+'.tsv'
with open(outPrf, 'w') as out: 
    i = 0
    while i < 100000:
        gene = all_gene_ids[np.random.randint(0, len(all_gene_ids))]
        
        seq = sequences[gene]
        if len(seq) <=  2 * window:  
            continue
        all_locs = []
        if gene in id2prfloc.keys():
             
            prfs = id2prfloc[gene]
       
            for p in prfs:
                all_locs.append(int(p))
         
        rand = np.random.randint(window, len(seq) - window)
        rand = rand - (rand % 3)
        while rand in all_locs:
            rand = np.random.randint(window, len(seq) - window)
            rand = rand - (rand % 3)
        
        
        before = seq[(rand - window):rand]
        after = seq[(rand): (rand + window)]

        cai_before, cost_before = calc_cai_cost(before.upper())
        cai_after,cost_after  = calc_cai_cost(after.upper())
        
        out.write(gene + '\t' + str(rand) + '\t' + str(len(seq)) + '\t' + str(cai_before) + '\t' + str(cai_after) + '\t'+ str(cai_after - cai_before) + '\t' + str(cost_before) +'\t' + str(cost_after) + '\t'+ str(round((cost_after - cost_before),4))+'\n') 
        i += 1 
        
            
        
        
        
        
            
        
        
        