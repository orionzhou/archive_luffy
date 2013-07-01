#!/bin/bash

cd $data/genome/Mtruncatula_4.0/raw

seqconcat.pl -i JCVI.Medtr.v4.unplaced.fasta -o JCVI.Medtr.v4.unplaced.concat.fasta -p JCVI.Medtr.v4.unplaced.concat.tbl -i chrU -g 1000

cat JCVI.Medtr.v4.chrs.fasta JCVI.Medtr.v4.unplaced.concat.fasta > ../11_genome.fa

gff_convert_loc.pl -i 11.gff -p JCVI.Medtr.v4.unplaced.concat.tbl -o 15_global_loc.gff
gff_to_gtb.pl -i 15_global_loc.gff -o ../21_gene.gtb -s ../11_genome.fa

cd ..
gtb_conv.pl -i 21_gene.gtb -o 21_gene.gff
gtb_conv.pl -i 21_gene.gtb -o 21_gene.fas -f seq2 -s ../11_genome.fa
gtb_conv.pl -i 21_gene.gtb -o 21_gene.tbl -f tbl