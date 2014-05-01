#!/bin/bash

cd $data/genome/HM056 #HM340

cd augustus
gff_augustus.pl -i 01.gff -o - | gff2gtb.pl -i - -o 02.gtb
gtb_augustus.pl -i 02.gtb -o 03.gtb
awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR!=1) {$17="gene";} print}' 03.gtb > ../40_gene.gtb
cd ..

awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR!=1) {$18=$17; $17="CRP"; print}}' spada.crp.gtb | cut -f1-18 > 43_crp.gtb

cat 4[03]*.gtb > 49.gtb
gtbdedup.pl -i 49.gtb -o 51.gtb


gtb2gff.pl -i 51.gtb -o 51.gff
gtb2tbl.pl -i 51.gtb -o 51.tbl
gtb2bigbed.pl -i 51.gtb -s 15.sizes -o 51.bb
#gtb2fas.pl -i 51.gtb -o 51.fas -s 11_genome.fa
