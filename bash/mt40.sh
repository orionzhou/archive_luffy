#!/bin/bash

cd $data/genome/Mtruncatula_4.0/raw

ln -sf JCVI.Medtr.v4.20130313.fasta 01_all.fas
grep '>' 01_all.fas | grep -o 'chr[0-9]\+' > 02.chrs.txt
grep '>' 01_all.fas | grep -o 'scaffold[0-9]\+' > 02.scaffolds.txt
seqret.pl -d 01_all.fas -b 02.chrs.txt -o 03.chrs.fas
seqret.pl -d 01_all.fas -b 02.scaffolds.txt -o 03.scaffolds.fas
seqconcat.pl -i 03.scaffolds.fas -o 03.scaffolds.concat.fas -p 03.scaffolds.concat.tbl -d chrU -g 1000

cat 03.chrs.fas 03.scaffolds.concat.fas > ../11_genome.fas

gff_jcvi.pl -i 11.gff -o 12_jcvi_fixed.gff
gff_convert_loc.pl -i 12_jcvi_fixed.gff -p 03.scaffolds.concat.tbl -o 15_global_loc.gff
gff2gtb.pl -i 15_global_loc.gff -o 15_global_loc.gtb
gtbcheckphase.pl -i 15_global_loc.gtb -s ../11_genome.fa -o 17_phase_fixed.gtb
gtbcatte.pl -i 17_phase_fixed.gtb -o 21.gtb
gtbdedup.pl -i 21.gtb -o 25_dedup.gtb
gtblongest.pl -i 25_dedup.gtb -o 26_longest.gtb

cd ..
ln -sf raw/26_longest.gtb 21.gtb
