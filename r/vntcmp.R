require(plyr)
require(VennDiagram)

org.q = "HM056"
dir.q = file.path('/home/youngn/zhoup/Data/genome', org.q)

fs.q = file.path(dir.q, '15_seqlen.tbl')
fg.q = file.path(dir.q, '16_gaploc.tbl')
ts.q = read.table(fs.q, sep='\t', header=T, as.is=T)
tg.q = read.table(fg.q, sep='\t', header=T, as.is=T)

org.t = "HM101"
dir.t = file.path('/home/youngn/zhoup/Data/genome', org.t)

fs.t = file.path(dir.t, '15_seqlen.tbl')
fg.t = file.path(dir.t, '16_gaploc.tbl')
ts.t = read.table(fs.t, sep='\t', header=T, as.is=T)
tg.t = read.table(fg.t, sep='\t', header=T, as.is=T)


# SNP density
fw = sprintf('/home/youngn/zhoup/Data/misc3/%s_%s/41_novelseq/nov2/35.gal', org.q, org.t)
tw = read.table(fw, sep='\t', header=T, as.is=T)[,1:17]
tw = cbind(tw, alnlen=tw$match+tw$misMatch)
tw = tw[tw$alnlen > 0,]

sum(tw$alnlen)
sum(tw$misMatch)
sum(tw$misMatch) / sum(tw$alnlen)

intvs = c(0,100,1000,10000,100000,Inf)
labels = c("1-100", "100-1000", "1k-10k", "10k-100k", "100k+")

tw.2 = cbind(tw, sub=tw$misMatch/tw$alnlen, intv=cut(tw$alnlen, breaks=intvs, labels=labels))
tws = ddply(tw.2, .(intv), summarise, cnt=length(sub), misMatch=sum(misMatch), alnlen=sum(alnlen), sub.mean=sum(misMatch)/(sum(match)+sum(misMatch)), sub.median=median(sub))

# compare SNPs called from two sources
chrs = seq(8)
names(chrs) = sprintf("chr%s", seq(8))

fd = sprintf('/home/youngn/zhoup/Data/misc3/%s_%s/41_novelseq/nov2/35.snp', org.q, org.t)
td = read.table(fd, sep='\t', header=T, as.is=T)[,1:2]
colnames(td) = c('chr', 'pos')
td.1 = td[td$chr %in% names(chrs),]
possd = chrs[td.1$chr]*1000000000+td.1$pos

fb = sprintf('/home/youngn/zhoup/Data/misc3/hapmap_mt40/30_vnt/%s.snp', org.q)
tb = read.table(fb, sep='\t', header=T, as.is=T)[,1:2]
colnames(tb) = c('chr', 'pos')
tb.1 = tb[tb$chr %in% names(chrs),]
possb = chrs[tb.1$chr]*1000000000+tb.1$pos

area1 = length(possd)
area2 = length(possb)
areac = sum(possd %in% possb)
venn.plot <- draw.pairwise.venn( area1, area2, areac, 
category = c("de novo", "read mapping"),
fill = c("blue", "red"), lty = "blank", cex = 2, 
cat.cex = 2, cat.pos = c(285, 105), cat.dist = 0.09, cat.just = list(c(-1, -1), c(1, 1)),
ext.pos = 30, ext.dist = -0.05, ext.length = 0.85, ext.line.lwd = 2, ext.line.lty = "dashed")
tiff(sprintf('/home/youngn/zhoup/Data/misc3/%s_%s/snpcmp.tiff', org.q, org.t))
grid.draw(venn.plot)
dev.off()


