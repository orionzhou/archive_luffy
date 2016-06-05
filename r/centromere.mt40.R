require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(gtable)
require(grid)
require(RColorBrewer)
source('Location.R')
source('comp.fun.R')

dirw = file.path(Sys.getenv("misc2"), "centromere")

tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2, stringsAsFactors=F)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3, stringsAsFactors=F)
tp = cbind(tp, len = tp$end-tp$beg+1)

### find largest gaps - not working
xf <- function(df) {
  df = cbind(df, len = df[,'end']-df[,'beg']+1)
  idxs = which(df[,'len']==max(df[,'len']))
  data.frame(chr=df[idxs,'chr'], beg=df[idxs,'beg'], end=df[idxs,'end'], len=df[idxs,'len'], stringsAsFactors=F)
}
tl = ddply(tp, .(chr), xf)
tl = tl[tl$chr != 'chrU',]
fl = file.path(dirw, "31.lgap.mt40.tbl")
write.table(tl, fl, sep = "\t", row.names = F, col.names = T, quote = F)

### using Mt3.5 centromeres
fi = file.path(Sys.getenv("genome"), "Mtruncatula_3.5/sequence/07_cen.tbl")
ti = read.table(fi, header = F, sep = "\t", as.is = T)
dc1 = within(ti, {V3=V2-1; V2 = V2-1000;})
dc2 = within(ti, {V2=V3+1; V3 = V3+1000;})
dm = rbind(dc1, dc2)
fm = file.path(dirw, "61.cen.mt35.bed")
write.table(dm, fm, sep = "\t", row.names = F, col.names = F, quote = F)
# seqret.pl -d $genome/Mtruncatula_3.5/11_genome.fa -b 61.cen.mt35.bed -o 61.fas
# blat $genome/HM101/11_genome.fas 61.fas 62.psl
# psl2gal.pl -i 62.psl -o 62.gal

fi = file.path(dirw, "62.gal")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
qChr = sapply(strsplit(ti$qId, "-"), "[", 1) 
tm = ti[ti$mat > 950 & ti$tId == qChr,]

### prepare sliding window of gap stat
x = tt$end
names(x) = tt$chr
gr = tileGenome(x, tilewidth = 100000, cut.last.tile.in.chrom = T)

tw = data.frame(chr = seqnames(gr), beg = start(gr), end = end(gr), 
  len = width(gr), stringsAsFactors = F)
bp_gap = intersect_basepair(gr, grp)
gap = bp_gap / tw$len

cols = brewer.pal(9, "Greys")

vmin = min(gap); vmax = max(gap)
breaks = seq(vmin, vmax, length.out = length(cols)+1)
gapcol = cut(gap, breaks, include.lowest = T)
tw = cbind(tw, gap=gap, gap=gapcol)


chrn = 1:nrow(tt)
names(chrn) = rev(tt$chr)

fi = file.path(dirw, "11.gal")
ti = read.table(fi, header = T, sep = "\t", as.is = T)

fc = file.path(dirw, "cen.mt40.tsv")
tc = read.table(fc, header = F, sep = "\t", as.is = T)
colnames(tc) = c('chr','beg','end')
tc = cbind(tc, pos = (tc$end+tc$beg)/2)

tt2 = cbind(tt, chrn = chrn[as.character(tt$chr)])
tw2 = cbind(tw, chrn = chrn[as.character(tw$chr)])
ti2 = cbind(ti, chrn = chrn[as.character(ti$tId)])
tm2 = cbind(tm, chrn = chrn[as.character(tm$tId)])
tc2 = cbind(tc, chrn = chrn[as.character(tc$chr)])

pc <- ggplot() +
  geom_rect(data = tw2, aes(xmin=beg, xmax=end, ymin=chrn-0.25, ymax=chrn+0.25, fill=gapcol), linetype=0) +
  geom_rect(data = tt2, aes(xmin=beg, xmax=end, ymin=chrn-0.25, ymax=chrn+0.25), fill=NA, color='deeppink4') +
  geom_point(data = ti2, aes(x=tBeg, y=chrn+0.4, col=qId, shape=qId), size=2) +
  geom_point(data = tm2, aes(x=tBeg, y=chrn+0.4, col='Mt3.5 Cen', shape='Mt3.5 Cen'), size=2) +
  geom_point(data = tc2, aes(x=pos, y=chrn-0.35, col='Predicted Mt4.0 Cen', shape='Predicted Mt4.0 Cen'), size=2) +
  theme_bw() + 
  scale_x_continuous(name = '', expand = c(0.01, 0)) + 
  scale_y_continuous(name = '', expand = c(0, 0), breaks=chrn, labels=names(chrn), limits=c(0.5,9.7)) +
  scale_fill_manual(name='', breaks=breaks, values=cols, guide=guide_legend(label.position='none')) +
  scale_color_manual(name='', values=c('red','blue','green','purple','black')) +
  scale_shape_manual(name='', values=c(6,6,6,6,17)) +
  theme(legend.position = c(0.9,0.45), legend.direction = "vertical", legend.title = element_text(size = 8), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size=8), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "cm")) +
#  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))

fo = file.path(dirw, "91.pdf")
ggsave(pc, filename = fo, width = 10, height = 6)


## look for centromeric repeats in HM340.PB assembly
tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2, stringsAsFactors=F)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3, stringsAsFactors=F)

fi = file.path(dirw, "51.HM340.gal")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti = ti[ti$score > 100 & ti$qId == 'MtR3',c(2:11,13:14,18:19)]

