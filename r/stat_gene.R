dir = file.path(DIR_Misc1, "stats_gene")
f_gc = file.path(dir, "01_gc.tbl")
t_gc = read.table(f_gc, header=T, as.is=T, sep="\t")
f_id = file.path(DIR_Misc1, "seq06/opt12/02_ids.tbl")
t_id = read.table(f_id, header=T, as.is=T, sep="\t")
t = merge(t_gc, t_id, by='id')
t.2 = t[,c('id','length.x','pct_gc')]
colnames(t.2)[2]='length'
write.table(t.2, file.path(dir, "02_gc_opt12.tbl"), sep="\t", quote=F, col.names=T, row.names=F)
