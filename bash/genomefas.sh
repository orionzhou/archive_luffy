#!/bin/bash

seqcheck.pl -i ***.fas -o 11_genome.fa

seqlen.pl -i 11_genome.fa -o 15.sizes
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, 0, $2}' 15.sizes > 15.bed

seqgap.pl -i 11_genome.fa -o 16_gap.bed -m 10
#awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR != 1) {print $1, $2-1, $3}}' \
#  16_gap.tbl > 16_gap.bed
bedToBigBed 16_gap.bed 15.sizes 16_gap.bb


