## Goals
1. Verify CAI I programmed with the CAI claculator internet tool found at https://www.biologicscorp.com/tools/CAICalculator/ -scores were identical so that is good news. 
2. Verify relative adaptation scores for S.cerevesia from webtool

## Steps
1. yeast_rel_adapt.txt contains the table of rel adapt from the website
2. Use split_rel_adapt.py to change yeast_rel_adapt.txt into a more usable format
3. Final file is yeast_rel_adapt_adj.txt
4. Use verify_yeast_colab.ipynb to calculate relative adaption scores in yeast by hand. Uses all single exon protein coding genes in yeast
5. Results comparing out calculation with the values from the webtool are found in compare_rel_adapt_scores.txt - overall very similar values
