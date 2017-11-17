require(dplyr)

dirw = '/home/springer/zhoux379/data/genome/Zmays_v4/v37'

fi = file.path(dirw, "t2.gtb")
ti = read.table(fi, sep = "\t", as.is = T, header = T)
colnames(ti)[2] = 'gid'

grp = dplyr::group_by(ti, gid)
tg = dplyr::summarise(grp, chr = chr[1], beg = min(beg), end = max(end), srd = srd[1])

to = tg[,c(2:4,1)]
to$beg = to$beg - 1
to = to[order(to$chr, to$beg),]

fo = file.path(dirw, "gene.bed")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
