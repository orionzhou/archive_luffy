require(dplyr)
require(GenomicRanges)
require(tidyr)

diri = file.path(Sys.getenv("misc2"), "briggs")
dirw = file.path(Sys.getenv("misc2"), "erika")

fx = file.path(diri, '35.long.tsv')
tx = read.table(fx, header = T, sep = "\t", as.is = T)

### TF expression
fi = file.path(dirw, 'tf_v4.tsv')
ti = read.table(fi, header = F, sep = "\t", as.is = T)

gids = ti$V2
gids = unique(gids[gids != ''])

tx2 = tx[tx$gid %in% gids,]

tx3 = cbind(tx2, treat = sprintf("%s|%s", tx2$Genotype, tx2$Tissue))
to1 = spread(tx3[,c('gid','treat','fpm')], treat, fpm)
to2 = spread(tx3[,c('gid','treat','fpkm')], treat, fpkm)
fo = file.path(dirw, "tf_v4_fpm.tsv")
write.table(to1, fo, sep = "\t", row.names = F, col.names = T, quote = F)

fo = file.path(dirw, "tf_v4_fpkm.tsv")
write.table(to2, fo, sep = "\t", row.names = F, col.names = T, quote = F)


### promoter expression
fi = file.path(dirw, '11.promoter.tsv')
ti = read.table(fi, header = F, sep = "\t", as.is = T)

gids = ti$V2
gids = unique(gids[gids != ''])
nrow(ti)
length(gids)

tx2 = tx[tx$gid %in% gids,]

tx3 = cbind(tx2, treat = sprintf("%s|%s", tx2$Genotype, tx2$Tissue))
to1 = spread(tx3[,c('gid','treat','fpm')], treat, fpm)
to2 = spread(tx3[,c('gid','treat','fpkm')], treat, fpkm)
fo = file.path(dirw, "12.promoter.fpm.tsv")
write.table(to1, fo, sep = "\t", row.names = F, col.names = T, quote = F)

fo = file.path(dirw, "12.promoter.fpkm.tsv")
write.table(to2, fo, sep = "\t", row.names = F, col.names = T, quote = F)
