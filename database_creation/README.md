## goal: 
create a subset of Saccharomyces cerevisiae S288C (baker's yeast) genes that are single exon and protein coding
## steps
1. Get gencode fasta and gff file from NCBI https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2 
2. Use remove_multis.py script to remove multiexon genes from the gff 
3. Use *grep gene_biotype=protein_coding yeast_gencode_single_exon.gff > out.gff* to keep only protein coding genes
4. Extract these gene sequences from the fasta file with the bedtool command: *bedtools getfasta -fi yeast.fasta -bed  yeast_single_exon_coding.gff -s > out.fa*  
