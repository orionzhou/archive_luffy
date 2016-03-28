require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(gtable)
require(grid)
require(gridExtra)
require(RColorBrewer)
#require(Gviz)
source('Location.R')
source('comp.win.fun.R')

dirw = file.path(Sys.getenv("misc3"), "comp.stat")

##### read in data tracks
tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

### gene fam stats
tg = read.table(tcfg$gene, sep = "\t", header = F, as.is = T)[,1:7]
colnames(tg) = c("chr", 'beg', 'end', 'srd', 'id', 'type', 'cat')
tg = tg[tg$type == 'cds',]

### covered bases
grl = list()
for (qname in qnames_15) {
  dirc = sprintf("%s/%s_%s", Sys.getenv("misc3"), qname, tname) 
  fa = file.path(dirc, '23_blat/31.9/gax')
  ta = read.table(fa, header = T, sep = "\t", as.is = T)
  colnames(ta) = c('tchr','tbeg','tend','tsrd','qchr','qbeg','qend','qsrd','cid','lev')
  gra = with(ta[ta$lev<=20,], GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend)))
  grl[[qname]] = reduce(gra)
}

### ingroup coverage & theta-pi
fc = file.path(Sys.getenv("misc3"), "comp.vnt", "81.cvg.tbl")
tc = read.table(fc, header = F, sep = "\t", as.is = T)
grc = with(tc, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

fn = file.path(Sys.getenv("misc3"), "comp.vnt", "25.stat.tbl")
tn = read.table(fn, header = T, sep = "\t", as.is = T)
grn = with(tn, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))

fv = file.path(Sys.getenv("misc3"), "comp.vnt", "52.stat.tbl")
tv = read.table(fv, header = T, sep = "\t", as.is = T)
tv = cbind(tv, size = (tv$rsize + tv$asize - 2))
tvs = tv[tv$size < 50,]
tvl = tv[tv$size >= 50,]
grvs = with(tvs, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))
grvl = with(tvl, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))

##### create sliding window table using 1Mb sliding windows
source('comp.win.fun.R')
chr = 'chr2'
beg = 7000001
title = sprintf("%s:%02dMb", chr, as.integer(beg/1000000))

res = prepare_data(chr, beg, grt, grp, tg, grl, grc, grn, grvs, grvl)
splots = sub_plots(chr, beg, res$tw, res$dg, res$dy, res$ds)

### multi-panel plot
gp1 = splots$cvg + theme(axis.text.y = element_blank())
gt1 = ggplotGrob(gp1)
#gt1 = gtable_add_cols(gt1, unit(0, "mm"))
gs = list(gt1)
heis = c(15)
for (key in names(splots$splots)) {
	gp2 = splots$splots[[key]] + theme(axis.text.y = element_blank())
	gt2 = ggplotGrob(gp2)
	gt2$widths = gt1$widths
	gs = c(gs, list(gt2))
	heis = c(heis, 1)
}
for (key in names(splots$gplots)) {
	gp2 = splots$gplots[[key]] + theme(axis.text.y = element_blank())
	gt2 = ggplotGrob(gp2)
	gt2$widths = gt1$widths
	gs = c(gs, list(gt2))
	heis = c(heis, 1)
}
g <- gtable_matrix(name = 'demo', grobs = matrix(gs, nrow = length(gs)), widths = 1, heights = heis)
pp <- gtable_matrix(name = 'demo', grobs = matrix(list(g, g), nrow = 1), widths = c(1, 1), heights = 1)

bog <- rectGrob(gp = gpar(col='black', fill=NA, lwd=2))
bo <- gtable_matrix(name = 'demo', grobs = matrix(list(bog, bog), nrow = 1), widths = c(1, 1), heights = 1)

fo = sprintf("%s/33.wins/%s.pdf", dirw, title)
pdf(file = fo, width = 6, height = 6, bg = 'transparent')
grid.newpage()
grid.draw(pp)
grid.draw(bo)
dev.off()

##### create multiple sliding-window plot
chrs = c('chr2', 'chr2', 'chr4', 'chr5', 'chr7')
begs = c(16, 30, 5, 6, 28)

