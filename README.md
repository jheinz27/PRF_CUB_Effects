# README
## Calculate CAI and Cost in Before and After Windows
**calc_cai_before_after_prf.py -w (window size)** - will calculate for one window

**run_all_windows.sh** -  is a bash script to run multiple windows and write them to directory "all". Use "./run_all_windows.sh" to run 

**visualize_all_final.ipynb** - is a jupyter notebook to graph the data

## Predictions

**predict_prfs.py -fi (fasta file ) -b (bp of buffer for making ranges) -w (window size bp) -o (output.tsv)** - Return regions of high CAI scores aka a prediction of a PRF site

**calc_accuracy.py < prediction_output.txt** - Calculate the accuracies of the predicted regions and how many known PRFs were found

## mRNAs of Similar Lengths 
**comp_same_lens.ipynb** - used to compare CAI scores of genes of similar lengths with and without prf sites

## Gene Expression

**gene_exp.ipynb** - investigate how CAI relates to gene expression  

## CUB After PRF Site
**after_prf.ipynb** - Graphs codon by codon frequency shortly before and after a prf site
