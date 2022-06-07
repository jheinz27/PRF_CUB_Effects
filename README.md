# Summer22
## stuff I've done
1. filtered multiexon genes out of yeast (s288) gencode annotation, only got rid of ~400 genes (can look into these later but want everything as simple as possible for now)
2. used bedtools to get sequences of genes out. 
   bedtools getfasta -fi yeast.fasta -bed  yeast_single_exon.gff -s > out.fa 
3. used split_rel_adapt.py to convert the relative adaptiveness table into a more usable format

## future goals: 
1. validate relative adaptive values by going though all gene seqs - bc who trusts the internet **DONE**
2. calc CAI over sliding windows for all yeast genes- try to detect pattern- Might be sus due to gencode annotaton. can use refseq as it's more conservative, but that seems a world problem. 
3. graph all the CAI's of the genes - need to figure out how to plot these on one graph bc of varying gene lengths. 
4. access prfdb for yeast
5. code in frameshift at the -1 prf and calculate all CAI again - expect a decline around the prf loction
