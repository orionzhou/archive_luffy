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

### find largest gaps
xf <- function(df) {
  idxs = which(df[,'len']==max(df[,'len']))
  data.frame(chr=df[idxs,'chr'], beg=df[idxs,'beg'], end=df[idxs,'end'], len=df[idxs,'len'], stringsAsFactors=F)
}
tc = ddply(tp, .(chr), xf)
tc = tc[tc$chr != 'chrU',]
fc = file.path(dirw, "31.cent.mt40.tbl")
write.table(tc, fc, sep = "\t", row.names = F, col.names = T, quote = F)

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

tt2 = cbind(tt, chrn = chrn[as.character(tt$chr)])
tw2 = cbind(tw, chrn = chrn[as.character(tw$chr)])
ti2 = cbind(ti, chrn = chrn[as.character(ti$tId)])
tc2 = cbind(tc, chrn = chrn[as.character(tc$chr)])

pc <- ggplot() +
  geom_rect(data = tw2, aes(xmin=beg, xmax=end, ymin=chrn-0.25, ymax=chrn+0.25, fill=gapcol), linetype=0) +
  geom_rect(data = tt2, aes(xmin=beg, xmax=end, ymin=chrn-0.25, ymax=chrn+0.25), fill=NA, color='deeppink4') +
  geom_point(data = ti2, aes(x=tBeg, y=chrn+0.45, col=qId), shape=6, size=2) +
  geom_point(data = tc2, aes(x=beg, y=chrn+0.45, col='gap'), shape=25, size=2) +
  theme_bw() + 
  scale_x_continuous(name = '', expand = c(0.01, 0)) + 
  scale_y_continuous(name = '', expand = c(0, 0), breaks=chrn, labels=names(chrn), limits=c(0.5,9.7)) +
  scale_fill_manual(name='', breaks=breaks, values=cols, guide=guide_legend(label.position='none')) +
  scale_color_brewer(name='', palette='Set1') +
  theme(legend.position = c(0.9,0.1), legend.direction = "vertical", legend.title = element_text(size = 8), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size=8), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "cm")) +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))

fo = file.path(dirw, "21.pdf")
ggsave(pc, filename = fo, width = 6, height = 4)


## look for centromeric repeats in HM340.PB assembly
tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2, stringsAsFactors=F)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3, stringsAsFactors=F)

fi = file.path(dirw, "51.HM340.gal")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti = ti[ti$score > 100 & ti$qId == 'MtR3',]

