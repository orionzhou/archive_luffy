#!/bin/bash

cd $data/genome/HM056
#cd $data/genome/HM340



awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR!=1) {$18=$17; $17="CRP"; print} else {print}}' spada.crp.gtb | cut -f1-18 > 43_crp.gtb

cat 4[0-3]*.gtb > 49.gtb
gtbdedup.pl -i 49.gtb -o 51.gtb

#awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR==1 || tolower($17) != "te") print}' 51.gtb > 55_noTE.gtb

gtb2gff.pl -i 51.gtb -o 51.gff
gtb2fas.pl -i 51.gtb -o 51.fas -s 11_genome.fa
gtb2bigbed.pl -i 51.gtb -s 15.sizes -o 51.bb
