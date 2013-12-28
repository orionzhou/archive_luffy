#!/bin/bash

cd $data/genome/Mtruncatula_4.0/raw

ln -sf JCVI.Medtr.v4.20130313.fasta 01_all.fa
grep '>' 01_all.fa | grep -o 'chr[0-9]\+' > 02.chrs.txt
grep '>' 01_all.fa | grep -o 'scaffold[0-9]\+' > 02.scaffolds.txt
seqextract.pl -o 03.chrs.fa -i 02.chrs.txt 01_all.fa
seqextract.pl -o 03.scaffolds.fa -i 02.scaffolds.txt 01_all.fa
seqconcat.pl -i 03.scaffolds.fa -o 03.scaffolds.concat.fa -p 03.scaffolds.concat.tbl -d chrU -g 1000

cat 03.chrs.fa 03.scaffolds.concat.fa > ../11_genome.fa

gff_jcvi.pl -i 11.gff -o 12_jcvi_fixed.gff
gff_convert_loc.pl -i 12_jcvi_fixed.gff -p 03.scaffolds.concat.tbl -o 15_global_loc.gff
gffcheckphase.pl -i 15_global_loc.gff -s ../11_genome.fa -o 17_phase_fixed.gff
gff2gtb.pl -i 17_phase_fixed.gff -o ../21_gene.gtb

cd ..
gtb2gff.pl -i 21_gene.gtb -o 21_gene.gff
gtb2fas.pl -i 21_gene.gtb -o 21_gene.fas -s 11_genome.fa
gtb2bigbed.pl -i 21_gene.gtb -o 21_gene.bb
gtb2tbl.pl -i 21_gene.gtb -o 21_gene.tbl

awk 'BEGIN {FS="\t"} {if(NR==1 || $15 == "gene") print}' 21_gene.gtb > 40_gene.gtb
awk 'BEGIN {FS="\t"} {if(NR==1 || $15 == "transposable_element_gene") print}' 21_gene.gtb > 41_te.gtb
awk 'BEGIN {FS="\t"} {if(tolower($18) ~ /nbs-lrr/) print}' 21_gene.gtb > 42_nbs.gtb

grep -i "" 21_gene.gtb > 42_nbs.gtb
