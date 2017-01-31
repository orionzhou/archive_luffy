require(dplyr)
require(plyr)
require(GenomicRanges)

dirw = file.path(Sys.getenv("genome"), "Zmays_v4/TF")

### 
fg = file.path(dirw, "../51.gtb")
tg = read.table(fg, header = T, sep = "\t", as.is = T)[,1:3]

fm = file.path(dirw, "../gene_mapping/maize.v3TOv4.geneIDhistory.txt")
tm = read.table(fm, header = F, sep = "\t", as.is = T)
colnames(tm) = c("ogid", 'ngid', "change", "method", "type")

fi = file.path(dirw, '01.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
colnames(ti) = c("gid", "tid", 'fam')
ti = ti[ti$gid != '',]
ti = ddply(ti, .(gid), summarise, fam = names(sort(table(fam), decreasing = T))[1])

tm2 = tm[tm$ogid %in% ti$gid,]
table(tm2$type)
tm3 = tm2[tm2$type == '1-to-1',]
table(tm3$change)
tm4 = tm3[tm3$ngid %in% tg$par,]
tm5 = merge(tm4, ti, by.x = 'ogid', by.y = 'gid')
stopifnot(nrow(tm4) == nrow(tm5))

fo = file.path(dirw, "11.TF.txt")
write(tm4$ngid, fo)


### find house keeping genes
dirw = file.path(Sys.getenv("genome"), "Zmays_v4/housekeeping")

fi = file.path(dirw, '01.tsv')
ti = read.table(fi, header = F, sep = "\t", as.is = T)
colnames(ti) = c("ogid", 'name')

ti2 = merge(ti, tm, by = 'ogid')
sum(ti2$ngid %in% tg$par)

fo = file.path(dirw, "11.tsv")
write.table(ti2, fo, sep = "\t", row.names = F, col.names = T, quote = F)