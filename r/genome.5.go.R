require(GenomicRanges)
require(dplyr)

dirg = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro'

fg = file.path(dirg, "../51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,1:2]
colnames(tg) = c("tid", "gid")

fm = file.path(dirg, "../57.longest.tsv")
tm = read.table(fm, sep = "\t", header = F, as.is = T)
colnames(tm) = c("gid", 'tid')

ft = file.path(dirg, "07.tsv")
tt = read.table(ft, sep = "\t", header = F, as.is = T)
colnames(tt) = c("tid", "gos")

tt2 = merge(tg, tt, by = 'tid', all.x = T)
stopifnot(nrow(tg) == nrow(tt2))
tt2 = tt2[order(tt2$gid, tt2$tid),]

fo1 = file.path(dirg, "10_mrna.tsv")
write.table(tt2[,c(1,3)], fo1, sep = "\t", row.names = F, col.names = F, quote = F, na = '')

tm2 = merge(tm, tt, by = 'tid', all.x = T)
stopifnot(nrow(tm2) == nrow(tm))
tm2 = tm2[order(tm2$gid, tm2$tid),]
fo2 = file.path(dirg, "11_gene.tsv")
write.table(tm2[,c(2,3)], fo2, sep = "\t", row.names = F, col.names = F, quote = F, na = '')
