require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(gtable)
#require(Gviz)
source('Location.R')
source('comp.fun.R')

dirw = file.path(Sys.getenv("misc3"), "comp.stat")

tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

##### create sliding window table
x = tt$end
names(x) = tt$chr
gr = tileGenome(x, tilewidth = 100000, cut.last.tile.in.chrom = T)

tw = data.frame(chr = seqnames(gr), beg = start(gr), end = end(gr), 
  len = width(gr), stringsAsFactors = F)
  
bp_gap = intersect_basepair(gr, grp)
bp_nogap = tw$len - bp_gap

to = cbind(tw, len_ng = bp_nogap)
fo = file.path(diro, "31.win.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### test correlation of pi ~ gene density
fw = file.path(dirw, "32.win.stat.tbl")
tw = read.table(fw, header = T, sep = "\t", as.is = T)
gr = with(tw, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

tg = read.table(tcfg$gene, sep = "\t", header = F, as.is = T)[,1:7]
colnames(tg) = c("chr", 'beg', 'end', 'srd', 'id', 'type', 'cat')
tg = tg[tg$type == 'cds',]
tgg = tg[tg$cat != 'TE',]
tgt = tg[tg$cat == 'TE',]
tgn = tg[tg$cat == 'NBS-LRR',]
tgc = tg[tg$cat == 'CRP',]
ggg = with(tgg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
ggt = with(tgt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
ggn = with(tgn, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
ggc = with(tgc, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

p_g = intersect_basepair(gr, reduce(ggg)) / tw$len_ng
p_t = intersect_basepair(gr, reduce(ggt)) / tw$len_ng
p_n = intersect_basepair(gr, reduce(ggn)) / tw$len_ng
p_c = intersect_basepair(gr, reduce(ggc)) / tw$len_ng
tgd = data.frame('gene' = p_g, 'te' = p_t, 'nbs' = p_n, 'crp' = p_c, as.is = T)

tcc = cbind(tw, p_g = p_g, p_t = p_t, p_n = p_n, p_c = p_c)
tcc = tcc[tcc$lenc > 10000,]
summary(lm(pi_snp ~ p_g + p_t + p_n + p_c, data = tcc))
summary(lm(pi_indel ~ p_g + p_t + p_n + p_c, data = tcc))
summary(lm(pi_sv ~ p_g + p_t + p_n + p_c, data = tcc))

stats = list()
for (cname in c('pi_snp', 'pi_indel', 'pi_sv')) {
  txts = c()
  for (gname in c('p_g', 'p_t', 'p_n', 'p_c')) {
    cc = cor(tcc[,cname], tcc[,gname])
    pv = cor.test(tcc[,cname], tcc[,gname])$p.value
    txt = sprintf("r = %.03f, p = %s", cc, prettyNum(pv, digits = 3))
    txts = c(txts, txt)
  }
  stats[[cname]] = matrix(txts, nrow = 1, 
    dimnames = list(NULL, c("Non-TE genes", "TEs", "NBS-LRR", "CRP")))
}

ds = t(do.call(rbind.data.frame, stats))
fo = file.path(diro, "35_pi_cc.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)


# test correlation of gene density with distance to centromere
ts = tw[tw$chr == 'chr5',]
to = transform(ts, pct_cds = bp_cds/bp_nogap, bp_dist = abs( (beg+end)/2 - 21450000))
fit = lm(pct_cds ~ bp_dist, data = to)
summary(fit)
plot(to$bp_dist, to$pct_cds)

# test non-NCR v.s. NCR
fb = file.path(Sys.getenv("genome"), tname, "51.gtb")
tb = read.table(fb, sep = "\t", header = T, as.is = T)[,c(1,17)]
tg2 = merge(tg, tb, by = 'id')

tgc1 = tg2[tg2$cat == 'CRP' & tg2$cat3 <= 'CRP1030',]
tgc2 = tg2[tg2$cat == 'CRP' & tg2$cat3 > 'CRP1030' & tg2$cat3 < 'CRP1600',]
tgc3 = tg2[tg2$cat == 'CRP' & tg2$cat3 >= 'CRP1600',]
ggc1 = with(tgc1, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
ggc2 = with(tgc2, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
ggc3 = with(tgc3, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

p_c1 = intersect_basepair(gr, reduce(ggc1)) / bp_nogap
p_c2 = intersect_basepair(gr, reduce(ggc2)) / bp_nogap
p_c3 = intersect_basepair(gr, reduce(ggc3)) / bp_nogap

pc = p_c1 + p_c2 + p_c3
fit <- lm(tw$pi ~ tw$gen + tw$tre + tw$nbs + p_c2)
summary(fit)

#### global comparison plot
grl = list()
for (qname in qnames_15) {
  dirc = sprintf("%s/%s_%s", Sys.getenv("misc3"), qname, tname)
  
  fa = file.path(dirc, '23_blat/31.9/gax')
  ta = read.table(fa, header = T, sep = "\t", as.is = T)
  colnames(ta) = c('tchr','tbeg','tend','tsrd','qchr','qbeg','qend','qsrd','cid','lev')
  gra = with(ta[ta$lev<=20,], GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend)))
  grl[[qname]] = gra
}

x = tt$end
names(x) = tt$chr
gr = tileGenome(x, tilewidth = 1000000, cut.last.tile.in.chrom = T)

tw = data.frame(chr = seqnames(gr), beg = start(gr), end = end(gr), 
  len = width(gr), stringsAsFactors = F)
  
bp_gap = intersect_basepair(gr, grp)
bp_nogap = tw$len - bp_gap

tw = cbind(tw, len_ng = bp_nogap)

tw = tw[tw$chr != 'chrU' & tw$len_ng/tw$len > 0.5,]
gr = with(tw, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

goff = cumsum(tt$end + 3000000) + 1
goff = c(1, goff[1:(length(goff)-1)])
names(goff) = tt$chr
tx = cbind(tt, gpos = goff + (tt$beg + tt$end) / 2)

tw = cbind(tw, gbeg = tw$beg + goff[tw$chr] - 1, gend = tw$end + goff[tw$chr] - 1)

to = data.frame()
for (qname in qnames_15) {
  len_syn = intersect_basepair(gr, reduce(grl[[qname]]))
  pc = len_syn / tw$len_ng
  to = rbind(to, cbind(tw, org = qname, len_syn = len_syn, pc = pc))
}
to$org = factor(to$org, levels = rev(qnames_15))

breaks = seq(max(to$pc), max(to$pc), length.out=12)
to = cbind(to, pci = cut(to$pc, breaks, include.lowest = T))

cols = rev(c("#282b68", "#324387", "#385193", "#498bbd", "#71c5cd", "#81c185", "#afcf45", "#faed29", "#ea862d", "#db382b", "#bb242a"))
labs = sort(unique(to$pci))

pb <- ggplot(to) +
  geom_tile(aes(x = gbeg, y = org, fill = pci, height = 0.8)) +
  theme_bw() + 
  scale_x_continuous(name = '', expand = c(0, 0), breaks = tx$gpos, labels = tx$chr) + 
  scale_y_discrete(expand = c(0, 0), name = '') +
#  scale_fill_gradient(name = '% covered in synteny alignment', space = "Lab", low = 'firebrick1', high = 'dodgerblue') +
#  scale_fill_distiller(name = '% covered in synteny alignment', type = "seq", space = "Lab", direction = 1, palette = "Spectral") +
  scale_fill_manual(name = '% covered in synteny alignment', breaks = labs, labels = labs, values = cols, guide = guide_legend(nrow = 2, byrow = T)) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0, 1), legend.title = element_text(size = 8), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,1,0,0), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "blue", size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))

fo = file.path(dirw, "12.comp.pdf")
ggsave(pb, filename = fo, width = 8, height = 4)

##### sliding window analysis plot
fw = file.path(dirw, "32.win.stat.tbl")
tw = read.table(fw, header = T, sep = "\t", as.is = T)
gr = with(tw, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

chr = 'chr5'
idxs = which(tw$chr == chr)
to = tw[idxs,]
to$pi_snp[to$lenc < 5000] = NA
to$pi_indel[to$lenc < 5000] = NA
to$pi_sv[to$lenc < 5000] = NA
tps = tp[tp$chr == chr,]

gd = tgd[idxs,]

chr_title = sprintf("%s position /Mbp", chr)
labs = seq(0, floor(tt$end[tt$chr == chr]) / 1000000, by = 10)
pb <- ggplot(cbind(to, gd)) +
  theme_bw() + 
  scale_x_continuous(name = chr_title, expand = c(0, 0), breaks = labs*1000000, labels = labs) + 
  theme(plot.margin = unit(c(0,1,0,0), "lines")) +
  theme(axis.title.x = element_text(colour = "black", size = 9)) +
  theme(axis.text.x = element_text(colour = "blue", size = 8)) +
  theme(axis.title.y = element_text(colour = "blue", size = 9)) +
  theme(axis.text.y = element_text(colour = "grey", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))

p_ng <- pb +
  geom_rect(data = tps, aes(xmin = beg, xmax = end, ymin = 0, ymax = 1)) +
  scale_y_continuous(expand = c(0, 0), name = 'Gap', breaks = NULL) + 
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())
p_gd1 <- pb +
  geom_rect(aes(xmin = beg, xmax = end, ymin = 0, ymax = gene, fill = 'Gene')) +
  geom_rect(aes(xmin = beg, xmax = end, ymin = gene, ymax = gene+te, fill = 'TE')) +
  theme(legend.position = c(0, 1.1), legend.direction = "horizontal", legend.justification = c(0, 1), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "line")) +
  scale_y_continuous(expand = c(0, 0), name = 'Gene density')
p_gd2 <- pb +
  geom_rect(aes(xmin = beg, xmax = end, ymin = 0, ymax = nbs, fill = 'NBS-LRR')) +
  geom_rect(aes(xmin = beg, xmax = end, ymin = nbs, ymax = nbs+crp, fill = 'CRP')) +
  theme(legend.position = c(0, 1.1), legend.direction = "horizontal", legend.justification = c(0, 1), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "line")) +
  scale_y_continuous(expand = c(0, 0), name = 'Gene density')
p_cb <- pb +
  geom_rect(aes(xmin = beg, xmax = end, ymin = 0, ymax = lenc/len)) +
  scale_y_continuous(expand = c(0, 0), name = 'Covered bases')
p_ps <- pb +
  geom_rect(aes(xmin = beg, xmax = end, ymin = 0, ymax = pi_snp)) +
  scale_y_continuous(expand = c(0, 0), name = 'Pi (SNP)')
p_pi <- pb +
  geom_rect(aes(xmin = beg, xmax = end, ymin = 0, ymax = pi_indel)) +
  scale_y_continuous(expand = c(0, 0), name = 'Pi (InDel)')
p_pv <- pb +
  geom_rect(aes(xmin = beg, xmax = end, ymin = 0, ymax = pi_sv)) +
  scale_y_continuous(expand = c(0, 0), name = 'Pi (SV)')

gt_ng <- ggplot_gtable(ggplot_build(p_ng)) 
gt_gd1 <- ggplot_gtable(ggplot_build(p_gd1)) 
gt_gd2 <- ggplot_gtable(ggplot_build(p_gd2)) 
gt_cb <- ggplot_gtable(ggplot_build(p_cb)) 
gt_ps <- ggplot_gtable(ggplot_build(p_ps))
gt_pi <- ggplot_gtable(ggplot_build(p_pi))
gt_pv <- ggplot_gtable(ggplot_build(p_pv)) 

tracks = list(gt_ng, gt_gd1, gt_gd2, gt_cb, gt_ps, gt_pi, gt_pv)
trackheights = c(1, 5, 5, 5, 5, 5, 5)

pad = mean(trackheights) * 0.05
hts = as.vector(rbind(rep(pad, length(tracks)), trackheights, rep(pad, length(tracks))))
gt <- gtable(widths = unit(c(1, 10, 0.1), "null"), height = unit(c(hts, 2), "null"))

for (i in 1:length(tracks)) {
  gt1 = tracks[[i]]
  gt <- gtable_add_grob(gt, gt1[3, 2:3], i*3-1, 1)
  gt <- gtable_add_grob(gt, gt1[3, 4], i*3-1, 2)
  if(ncol(gt1) == 6) {
    gt <- gtable_add_grob(gt, gt1[3, 5], i*3-1, 3)
  }
}
gt <- gtable_add_grob(gt, gt1[4:5, 4], length(tracks)*3 + 1, 2)

fo = sprintf("%s/33.vnt.stat.%s.pdf", dirw, chr)
pdf(file = fo, width = 8, height = 8, bg = 'transparent')
grid.newpage()
grid.draw(gt)
dev.off()

## calculate SNP density tracks
x = tt$end
names(x) = tt$chr
gr = tileGenome(x, tilewidth = 100, cut.last.tile.in.chrom = T)

bp_gap = intersect_basepair(gr, grp)
tw = data.frame(chr = seqnames(gr), beg = start(gr), end = end(gr), 
  len = width(gr), stringsAsFactors = F)

org = "HM034"

dirv = sprintf("%s/%s_HM101/23_blat", Sys.getenv("misc3"), org)
fv = sprintf("%s/31.9/snp", dirv)
tv = read.table(fv, sep = '\t', header = F, as.is = T)
colnames(tv) = c("chr", "pos", "ref", "alt", "qid", "qpos", "cid", "lev")
tv = tv[tv$lev == 1, c(1:2,4)]
grv = with(tv, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))

fx = sprintf("%s/%s_HM101/23_blat/31.9/gax", Sys.getenv("misc3"), org)
tx = read.table(fx, sep = '\t', header = F, as.is = T)
gr_gax = with(tx[tx$V10==1,], GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

bp_gax = intersect_basepair(gr, gr_gax)
cnt_snp = intersect_count(gr, grv)

to = cbind(tw, len_ng = tw$len - bp_gap, len_gax = bp_gax, snpc = cnt_snp)
to = to[to$len_gax > to$len * 0.3,]
to = cbind(to, snpd = to$snpc / to$len_gax)

ton = within(to, { beg = beg - 1; rm(len, len_ng, len_gax, snpc) })
fo = sprintf("%s/31.9/pct.bg", dirv)
options(scipen = 999)
write.table(ton, file = fo, sep = "\t", row.names = F, col.names = F, quote = F, na = '')
options(scipen = 0)

