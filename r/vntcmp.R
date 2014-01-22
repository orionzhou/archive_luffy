require(plyr)
require(VennDiagram)

qorg = "HM340"
qdir = file.path('/home/youngn/zhoup/Data/genome', qorg)
qfs = file.path(qdir, '15.sizes')
qfg = file.path(qdir, '16_gap.tbl')
qts = read.table(qfs, sep='\t', header=T, as.is=T)
qtg = read.table(qfg, sep='\t', header=T, as.is=T)

torg = "HM101"
tdir = file.path('/home/youngn/zhoup/Data/genome', torg)
tfs = file.path(tdir, '15.sizes')
tfg = file.path(tdir, '16_gap.tbl')
tts = read.table(tfs, sep='\t', header=T, as.is=T)
ttg = read.table(tfg, sep='\t', header=T, as.is=T)


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

dir = sprintf("%s/%s_%s/23_blat", DIR_Misc3, qorg, torg)
fd = file.path(dir, "27.snp")
td = read.table(fd, sep='\t', header=T, as.is=T)[,3:4]
colnames(td) = c('chr', 'pos')
td.1 = td[td$chr %in% names(chrs),]
possd = chrs[td.1$chr]*1000000000+td.1$pos

fb = sprintf('/home/youngn/zhoup/Data/misc3/hapmap_mt40/30_vnt/%s.snp', qorg)
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
tiff(sprintf('/home/youngn/zhoup/Data/misc3/%s_%s/snpcmp.tiff', qorg, torg))
grid.draw(venn.plot)
dev.off()


