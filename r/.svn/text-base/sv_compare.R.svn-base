a1 = readAssembly("mt_35")
a = drop.levels(subset(a1[a1$chr == "chr5",], type!='centromere'))
c = drop.levels(subset(a1[a1$chr == "chr5",], type=='centromere'))

acc1 = "HM010"
acc2 = "HM004"
accR = "HM101"
chr = "chr5"

#get coverage ratio analysis results
dir = file.path(DIR_Misc2, "cnv")
c1 = read.table(file.path(dir, "31.tbl"), sep="\t", head=T, stringsAsFactors=F)
c2 = c1[,c('id','chr','beg','end','p12_j','sig12','c12')]
delc12 = c2[!is.na(c2$sig12) & c2$sig12 == 1 & !is.na(c2$c12) & c2$c12 == 1,]
delc21 = c2[!is.na(c2$sig12) & c2$sig12 == 1 & !is.na(c2$c12) & c2$c12 == 2,]

#get PEM deletion results
dir = file.path(DIR_Repo, "mt_35/40_sv/41_shared")
d1 = read.table(file.path(dir, "11.tbl"), sep="\t", as.is=T, header=T, quote="")
d2 = d1[,c('id_pindel', 'chr', 'beg', 'end', 'size_d', 'size_i')]
g1 = read.table(file.path(dir, "71_genotype.tbl"), sep="\t", as.is=T, header=T, quote="")

del_id_12 = g1$id[g1[,acc1]=='G' & g1[,acc2]=='A']
del_id_21 = g1$id[g1[,acc1]=='A' & g1[,acc2]=='G']
delp12 = d2[d2$id_pindel %in% del_id_12,]
delp21 = d2[d2$id_pindel %in% del_id_21,]

#compare the two datasets
dir = file.path(DIR_Misc2, "cnv/compare")
write.table(delc12[,c('chr', 'beg', 'end', 'id')], file.path(dir, "delc12.bed"), sep="\t", quote=F, row.names=F, col.names=F)
write.table(delp12[,c('chr', 'beg', 'end', 'id_pindel')], file.path(dir, "delp12.bed"), sep="\t", quote=F, row.names=F, col.names=F)

write.table(delc21[,c('chr', 'beg', 'end', 'id')], file.path(dir, "delc21.bed"), sep="\t", quote=F, row.names=F, col.names=F)
write.table(delp21[,c('chr', 'beg', 'end', 'id_pindel')], file.path(dir, "delp21.bed"), sep="\t", quote=F, row.names=F, col.names=F)


