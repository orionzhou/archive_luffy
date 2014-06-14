
gtb2fas.pl -i 21.gtb -o 21.fas -d 11_genome.fas
# qsub interpro -N ipr.pfam -v FI=,FO=



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
