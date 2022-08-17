#!/bin/bash -xe 

for i in 15 30 60 90 120 150 300 450 600 750 900 1050 1200 1350 1500
do
	python calc_cai_before_after_prf.py -w $i 
done
mv all*.tsv all
