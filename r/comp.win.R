require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(gtable)
require(grid)
require(RColorBrewer)
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

##### create sliding window table using 1Mb sliding windows
x = tt$end
names(x) = tt$chr
gr = tileGenome(x, tilewidth = 1000000, cut.last.tile.in.chrom = T)

tw = data.frame(chr = seqnames(gr), beg = start(gr), end = end(gr), 
  len = width(gr), stringsAsFactors = F)
  
bp_gap = intersect_basepair(gr, grp)
bp_nogap = tw$len - bp_gap

tw = cbind(tw, len_ng = bp_nogap)
tw = tw[tw$chr != 'chrU' & tw$len_ng/tw$len > 0.5,]
fo = file.path(dirw, "31.win.tbl")
write.table(tw, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### gene family stats - # of genes per 100kb
fw = file.path(dirw, "31.win.tbl")
tw = read.table(fw, sep = "\t", header = T, as.is = T)
gr = with(tw, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

fg = file.path(tcfg$gdir, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1,3:5,16)]
colnames(tg) = c("id", "chr", 'beg', 'end', 'cat')

h = list('TE'='TE',
  'NBS' = c('TIR-NBS-LRR', 'CC-NBS-LRR'),
  'NCR' = 'NCR',
  'LRR_RLK' = 'LRR-RLK', 
  'F_box' = 'F-box'
)

to = tw[,1:3]
for (key in names(h)) {
  tgs = tg[tg$cat %in% h[[key]],]
  gs = with(tgs, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
  gd = intersect_count(gr, reduce(gs)) / tw$len_ng * 100000
  to = cbind(to, gd = gd)
  colnames(to)[length(to)] = key
}
fo = file.path(dirw, "32.win.gene.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### proportion covered based in synteny
grl = list()
for (qname in qnames_15) {
  dirc = sprintf("%s/%s_%s", Sys.getenv("misc3"), qname, tname) 
  fa = file.path(dirc, '23_blat/31.9/gax')
  ta = read.table(fa, header = T, sep = "\t", as.is = T)
  colnames(ta) = c('tchr','tbeg','tend','tsrd','qchr','qbeg','qend','qsrd','cid','lev')
  gra = with(ta[ta$lev<=20,], GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend)))
  grl[[qname]] = gra
}

to = tw[,1:3]
for (qname in qnames_15) {
  len_syn = intersect_basepair(gr, reduce(grl[[qname]]))
  pc = len_syn / tw$len_ng
  to = cbind(to, org = pc)
  colnames(to)[length(to)] = qname
}
fo = file.path(dirw, "32.win.syn.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### theta-pi, covered bases in 9/12 accessions
fc = file.path(Sys.getenv("misc3"), "comp.vnt", "81.cvg.tbl")
tc = read.table(fc, header = F, sep = "\t", as.is = T)
grc = with(tc, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
lenc = intersect_basepair(gr, grc)

### calc nuc-div for SNPs
fn = file.path(Sys.getenv("misc3"), "comp.vnt", "25.stat.tbl")
tn = read.table(fn, header = T, sep = "\t", as.is = T)
grn = with(tn, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))
nucdivs = intersect_score(gr, grn)
pi_snp = nucdivs / lenc

# (obsolete!) calc nuc-div for SNPs - the hard way
cl = makeCluster(detectCores())
cluster_fun <- function() {}
clusterCall(cl, cluster_fun)

get_nuc_div <- function(rw, fs) {
  chr = rw['chr']; beg = as.numeric(rw['beg']); end = as.numeric(rw['end'])
  cmd = sprintf("tabix %s %s:%d-%d | cut -f6 | paste -sd+ | bc", fs, chr, beg, end)
  ret = as.numeric(system(cmd, intern = T))
  ifelse(is.numeric(ret), ret, NA)
}
fs = file.path(dirw, '25.stat.tbl.gz')

ptm <- proc.time()
pis = parApply(cl, tw, 1, get_nuc_div, fs)
proc.time() - ptm

### calc nuc-div for indels/sv
fv = file.path(Sys.getenv("misc3"), "comp.vnt", "52.stat.tbl")
tv = read.table(fv, header = T, sep = "\t", as.is = T)
tv = cbind(tv, size = (tv$rsize + tv$asize - 2))
tvs = tv[tv$size < 50,]
tvl = tv[tv$size >= 50,]
grvs = with(tvs, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))
grvl = with(tvl, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))

nucdivs = intersect_score(gr, grvs)
pi_indel = nucdivs / tw$len_ng
nucdivs = intersect_score(gr, grvl)
pi_sv = nucdivs / tw$len_ng

pi_snp[is.infinite(pi_snp)] = NA
pi_indel[is.infinite(pi_indel)] = NA
pi_sv[is.infinite(pi_sv)] = NA
to = cbind(tw[,1:3], Coverage = lenc/tw$len_ng, Pi_SNP = pi_snp, Pi_InDel = pi_indel, Pi_SV = pi_sv)
fo = file.path(dirw, "32.win.stat.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

##### sliding window plot (1Mb)
fw = file.path(dirw, "31.win.tbl")
tw = read.table(fw, sep = "\t", header = T, as.is = T)

goff = cumsum(tt$end + 5000000) + 1
goff = c(1, goff[1:(length(goff)-1)])
names(goff) = tt$chr
tx = cbind(tt, gpos = goff + (tt$beg + tt$end) / 2, gbeg = goff + tt$beg - 1, gend = goff + tt$end - 1)
tx = tx[tx$chr != 'chrU',]

gbeg = tw$beg + goff[tw$chr] - 1
gend = tw$end + goff[tw$chr] - 1

fy = file.path(dirw, "32.win.syn.tbl")
ty = read.table(fy, sep = "\t", header = T, as.is = T)
ty = cbind(ty[,4:18], gbeg = gbeg, gend = gend)
tyl = reshape(ty, direction = 'long', varying = list(1:15), idvar = c('gbeg','gend'), timevar = 'org', v.names = 'syn', times = colnames(ty)[1:15])
tyl$org = factor(tyl$org, levels = rev(qnames_15))

cols1 = rev(c("#282b68", "#324387", "#385193", "#498bbd", "#71c5cd", "#81c185", "#afcf45", "#faed29", "#ea862d", "#db382b", "#bb242a"))
colsy = brewer.pal(11, "Spectral")

vmin = min(tyl$syn)
vmax = max(tyl$syn)
breaks = seq(vmin, vmax, length.out = length(colsy) + 1)
tyl = cbind(tyl, syn2 = cut(tyl$syn, breaks, include.lowest = T))

labs = sort(unique(tyl$syn2))
labg = rep('', length(colsy))
#labg[as.integer((1+length(colsy))/2)] = sprintf("%.03f - %.03f", vmin, vmax)
labg[length(colsy)] = sprintf("(%.03f-%.03f)", vmin, vmax)
ltitle = sprintf("proportion covered in synteny alignment (%.03f-%.03f)", vmin, vmax)

py <- ggplot(tyl) +
  geom_tile(aes(x = gbeg, y = org, fill = syn2, height = 0.8)) +
  geom_vline(data = tx, aes(xintercept = gbeg), size = 0.5) +
  geom_vline(data = tx, aes(xintercept = gend), size = 0.5) +
  theme_bw() + 
  scale_x_continuous(name = '', expand = c(0, 0), breaks = tx$gpos, labels = tx$chr) + 
  scale_y_discrete(expand = c(0, 0), name = '') +
#  scale_fill_gradient(name = '% covered in synteny alignment', space = "Lab", low = 'firebrick1', high = 'dodgerblue') +
#  scale_fill_distiller(name = '% covered in synteny alignment', type = "seq", space = "Lab", direction = 1, palette = "Spectral") +
  scale_fill_manual(name = ltitle, breaks = labs, labels = labg, values = colsy, guide = guide_legend(nrow = 1, byrow = T, label = T, label.position = 'top')) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0,0), legend.title = element_text(size = 8), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_blank(), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "line")) +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype = 0)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,1), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "blue", size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))

#fo = file.path(dirw, "12.comp.pdf")
#ggsave(py, filename = fo, width = 8, height = 4)

## gene fam density
fg = file.path(dirw, "32.win.gene.tbl")
tg = read.table(fg, sep = "\t", header = T, as.is = T)
tg = tg[,-1:-3]
fams = colnames(tg)
tg = cbind(tg, gbeg = gbeg, gend = gend)

lst = list()
for (fam in fams) {
colsg = brewer.pal(9, "Greens")

vmin = min(tg[,fam]); vmax = max(tg[,fam])
breaks = seq(vmin, vmax, length.out = length(colsg)+1)
tg2 = cbind(tg, gd = cut(tg[,fam], breaks, include.lowest = T))

labs = sort(unique(tg2$gd))
labg = rep('', length(colsg))
labg[length(colsg)] = sprintf("(%.03f-%.03f)", vmin, vmax)
ltitle = sprintf("number per 100kb (%.01f-%.01f)", vmin, vmax)

pg1 <- ggplot(tg2) +
  geom_tile(aes(x = gbeg, y = fam, fill = gd, height = 1)) +
  geom_segment(data=tx, aes(x=gbeg,xend=gend,y=0.5,yend=0.5), size = 0.5) +
  geom_segment(data=tx, aes(x=gbeg,xend=gend,y=1.5,yend=1.5), size = 0.5) +
  geom_vline(data = tx, aes(xintercept = gbeg), size = 0.5) +
  geom_vline(data = tx, aes(xintercept = gend), size = 0.5) +
  theme_bw() + 
  scale_x_continuous(name = '', expand = c(0, 0), breaks = tx$gpos, labels = tx$chr) + 
  scale_y_discrete(expand = c(0, 0), name = '') +
  scale_fill_manual(name = ltitle, breaks = labs, labels = labg, values = colsg, guide = guide_legend(nrow = 1, byrow = T, label = T, label.position = 'top')) +
  theme(legend.position = 'none', legend.direction = "horizontal", legend.justification = c(0,0), legend.title = element_text(size = 8), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_blank(), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "line")) + 
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype = 0)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))
lst[[fam]] = pg1
}

## theta-pi stats
fs = file.path(dirw, "32.win.stat.tbl")
ts = read.table(fs, sep = "\t", header = T, as.is = T)
ts = ts[,-1:-3]
ts = cbind(Gaps = 1-tw$len_ng/tw$len, ts)
keys = colnames(ts)
ts = cbind(ts, gbeg = gbeg, gend = gend)

stats = list()
for (key in keys) {
colst = brewer.pal(9, "Purples")
if(key == 'Coverage') {colst = brewer.pal(9, "Reds")}

vmin = min(ts[,key]); vmax = max(ts[,key])
breaks = seq(vmin, vmax, length.out = length(colst)+1)
ts2 = cbind(ts, val = cut(ts[,key], breaks, include.lowest = T))

labs = sort(unique(ts2$val))
labg = rep('', length(colst))
labg[length(colst)] = sprintf("(%.03f-%.03f)", vmin, vmax)
ltitle = sprintf("number per 100kb (%.01f-%.01f)", vmin, vmax)

ps <- ggplot(ts2) +
  geom_tile(aes(x = gbeg, y = key, fill = val, height = 1)) +
  geom_segment(data=tx, aes(x=gbeg,xend=gend,y=0.5,yend=0.5), size = 0.5) +
  geom_segment(data=tx, aes(x=gbeg,xend=gend,y=1.5,yend=1.5), size = 0.5) +
  geom_vline(data = tx, aes(xintercept = gbeg), size = 0.5) +
  geom_vline(data = tx, aes(xintercept = gend), size = 0.5) +
  theme_bw() + 
  scale_x_continuous(name = '', expand = c(0, 0), breaks = tx$gpos, labels = tx$chr) + 
  scale_y_discrete(expand = c(0, 0), name = '') +
  scale_fill_manual(name = ltitle, breaks = labs, labels = labg, values = colst, guide = guide_legend(nrow = 1, byrow = T, label = T, label.position = 'top')) +
  theme(legend.position = 'none', legend.direction = "horizontal", legend.justification = c(0,0), legend.title = element_text(size = 8), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_blank(), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "line")) +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype = 0)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))
stats[[key]] = ps
}

### multi-panel plot
g1 = ggplotGrob(py)
#g1 = gtable_add_cols(g1, unit(0, "mm"))
gs = list(g1)
heis = c(15)
for (key in keys) {
	g3 = ggplotGrob(stats[[key]])
	g3$widths = g1$widths
	gs = c(gs, list(g3))
	heis = c(heis, 1)
}
for (fam in fams) {
	g2 = ggplotGrob(lst[[fam]])
	g2$widths = g1$widths
	gs = c(gs, list(g2))
	heis = c(heis, 1)
}
g <- gtable_matrix(name = 'demo', grobs = matrix(gs, nrow = length(gs)), widths = 1, heights = heis)
#g$layout[grepl("guide", g$layout$name),c("t","b")] <- c(1,nrow(g))

fo = file.path(dirw, "12.comp.pdf")
pdf(file = fo, width = 8.5, height = 6.5, bg = 'transparent')
grid.newpage()
grid.draw(g)
dev.off()


##### test correlation of theta-pi and gene density (100kb sliding windows)
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

