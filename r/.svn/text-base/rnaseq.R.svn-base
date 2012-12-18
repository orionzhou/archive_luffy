dir = file.path(DIR_Misc2, "rnaseq")
fi = file.path(dir, "31_exp.tbl")
t1 = read.table(fi, sep="\t", header=T, as.is=T)[,c(1,5,2:4,6,7)]

exp_max = apply(t1[,-1], 1, max)
idx1 = which(exp_max >= 0.001)
t2 = t1[idx1,]

exp_max_rn = apply(t2[,c('nod', 'root')], 1, max)
exp_max_other = apply(t2[,c('flow','blade','pod','bud')], 1, max)
idx2 = which(exp_max_rn > 100*exp_max_other)
t3 = t2[idx2,]
write.table(t3, file.path(dir, "32_nod_root.tbl"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

