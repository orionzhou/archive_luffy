#!/bin/bash

cd $data/genome/pan3

seqlen.pl -i 11_genome.fa -o 15.sizes
seqgap.pl -i 11_genome.fa -o 16_gap.tbl -m 100
awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR>1) {$2=$2-1; print $1, $2, $3}}' 16_gap.tbl > 16_gap.bed
bedToBigBed 16_gap.bed 15.sizes 16_gap.bb

ln -sf ../Mtruncatula_4.0/40_gene.gtb 40_gene.gtb
ln -sf ../Mtruncatula_4.0/41_te.gtb 41_te.gtb
ln -sf ../Mtruncatula_4.0/42_nbs.gtb 42_nbs.gtb

awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR!=1) {$18=$17; $17="CRP"; print}}' spada.crp.gtb | cut -f1-18 > 43_crp.gtb

cat 4[0-3]*.gtb > 49.gtb
gtbdedup.pl -i 49.gtb -o 51.gtb

awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR==1 || tolower($17) != "te") print}' 51.gtb > 55_noTE.gtb

gtb2gff.pl -i 51.gtb -o 51.gff
gtb2fas.pl -i 51.gtb -o 51.fas -s 11_genome.fa
gtb2bigbed.pl -i 51.gtb -s 15.sizes -o 51.bb
