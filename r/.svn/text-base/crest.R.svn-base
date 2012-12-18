accs = get_mt_ids("acc26")

#filter Pindel output 
dir = file.path(DIR_Repo, "mt_35/40_sv/31_pindel")
f1 = file.path(dir, "11.tbl")
p1 = read.table(f1, sep="\t", header=T, as.is=T)
p2 = p1[p1$size_d >= 30 & p1$n_ind >=2 & p1$n_reads_uniq >= 5,]
write.table(p2, file.path(dir, "12_filtered.tbl"), sep="\t", quote=F, row.names=F, col.names=T)

#process and filter CREST outputs
dir = file.path(DIR_Repo, "mt_35/40_sv/33_crest")
fi = file.path(dir, "11_sum.tbl")
t01 = read.table(fi, header=T, sep="\t", as.is=T)
t01 = t01[t01$acc %in% accs,]
t02 = aggregate(t01$acc, by=list(factor(t01$chr_l), factor(t01$pos_l), factor(t01$strand_l), factor(t01$chr_r), factor(t01$pos_r), factor(t01$strand_r), factor(t01$type)), FUN=strconcat)

colnames(t02) = c("chr_l", "pos_l", "strand_l", "chr_r", "pos_r", "strand_r", "type", "acc")
t02$pos_l = as.numeric(as.character(t02$pos_l))
t02$pos_r = as.numeric(as.character(t02$pos_r))
t03 = cbind(id=1:nrow(t02), t02, len = t02$pos_r-t02$pos_l-1, n_acc = as.numeric(lapply(strsplit(t02$acc, " "), FUN=length)))
write.table(t03, file.path(dir, "21.tbl"), sep="\t", quote=F, row.names=F, col.names=T)

t04 = t03[t03$n_acc >= 2 & t03$type == 'DEL' & t03$len < 5000, ]
write.table(t04, file.path(dir, "22_filtered.tbl"), sep="\t", quote=F, row.names=F, col.names=T)

#get overlap btw Pindel & CREST
dir = file.path(DIR_Repo, "mt_35/40_sv/41_shared")

fp = file.path(dir, "../31_pindel/12_filtered.tbl")
p1 = read.table(fp, sep="\t", header=T, as.is=T)
fc = file.path(dir, "../33_crest/22_filtered.tbl")
c1 = read.table(fc, sep="\t", header=T, as.is=T)

c2 = cbind(c1, id_pindel=NA)
for (i in 1:nrow(c2)) {
  p2 = p1[abs(p1$beg-c2$pos_l[i]) <= 10 & abs(p1$end-c2$pos_r[i]) <= 10,]
#  c2$id_pindel[i] = length(p2$id)
  if(nrow(p2) >= 1) {c2$id_pindel[i] = p2$id[which.max(p2$n_ind)]}
}
colnames(c2)[1] = "id_crest"
colnames(p1)[1] = "id_pindel"
c3 = merge(p1, c2, by="id_pindel")

idxs = c()
for (id_pindel in unique(c3$id_pindel)) {
  c4 = c3[c3$id_pindel == id_pindel,]
  idx = which.max(c4$n_acc)
  idxs = c(idxs, row.names(c4)[idx])
}
c5 = c3[idxs,]
write.table(c5, file.path(dir, "11.tbl"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(c5[,1:4], file.path(dir, "12_simple.tbl"), sep="\t", quote=F, row.names=F, col.names=T)


