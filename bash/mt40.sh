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
gff2gtb.pl -i 15_global_loc.gff -o 15_global_loc.gtb
gtbcheckphase.pl -i 15_global_loc.gtb -s ../11_genome.fa -o 17_phase_fixed.gtb
gtbcatte.pl -i 17_phase_fixed.gtb -o 21.gtb
gtbdedup.pl -i 21.gtb -o 25_dedup.gtb
gtblongest.pl -i 25_dedup.gtb -o 26_longest.gtb

cd ..

awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR==1 || $17 == "gene") print}' 21.gtb > 40_gene.gtb
awk 'BEGIN {FS="\t"; OFS="\t"} {if($17 == "TE") print}' 21.gtb > 41_te.gtb

awk 'BEGIN {FS="\t"; OFS="\t"} {if(tolower($18) ~ /nbs-lrr/ || tolower($18) ~ /nb-arc/) {$17="NBS"; print}}' 21.gtb > 42_nbs.gtb
awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR!=1) {$18=$17; $17="CRP"; print}}' spada.crp.gtb | cut -f1-18 > 43_crp.gtb

cat 4[0-3]*.gtb > 49.gtb
gtbdedup.pl -i 49.gtb -o 51.gtb

awk 'BEGIN {FS="\t"; OFS="\t"} {if(NR==1 || tolower($17) != "te") print}' 51.gtb > 55_noTE.gtb

gtb2gff.pl -i 51.gtb -o 51.gff
gtb2tbl.pl -i 51.gtb -o 51.tbl
gtb2bigbed.pl -i 51.gtb -s 15.sizes -o 51.bb
#gtb2fas.pl -i 51.gtb -o 51.fas -s 11_genome.fa
