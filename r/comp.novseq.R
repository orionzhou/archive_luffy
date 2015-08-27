require(plyr)
require(dplyr)
require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(RColorBrewer)
source("Location.R")
source("comp.fun.R")

dirw = file.path(Sys.getenv('misc3'), 'comp.panseq')
diro = file.path(Sys.getenv('misc3'), 'comp.stat')

##### compile contaminated scaffold IDs
qname = "HM340"
for (qname in qnames) {
cfg = cfgs[[qname]]
ccfg = ccfgs[[qname]]
tl = read.table(cfg$size, header = F, sep = "\t", as.is = T)
colnames(tl) = c("chr", "size")

fn = file.path(ccfg$cdir, "../41_novseq/12.foreign.bed")
tn = read.table(fn, header = F, sep = "\t", as.is = T)
colnames(tn) = c("chr", "beg", 'end', 'type')
tn$beg = tn$beg + 1
tn = cbind(tn, flen = tn$end - tn$beg + 1)

tm = ddply(tn, .(chr), summarise, flen = sum(flen))
tm = merge(tm, tl, by = 'chr')
tm = cbind(tm, fpct = tm$flen/tm$size)
fids = tm$chr[tm$fpct >= 0.5]
x1 = sum(tn$flen[tn$chr %in% fids]) / sum(tn$flen)
x2 = sum(tn$flen[!tn$chr %in% fids])

tg = read.table(cfg$gene, header = F, sep = "\t", as.is = T)
tg = tg[tg$V6 == 'mrna', c(1,5)]
colnames(tg) = c('chr', 'gid')
y = sum(tg$chr %in% fids)
cat(sprintf("%s: %d scfs rmvd [%d genes], %.03f nov-seq rmvd, %d bp left\n", qname, length(fids), y, x1, x2))

fo = file.path(ccfg$cdir, "../41_novseq/15.foreign.scf.txt")
write(fids, fo)
}

##### enrichment analysis (by proportion of gene family being novel)
do = data.frame()
for (qname in qnames_all[1:12]) {
cfg = cfgs[[qname]]

tg = read.table(cfg$gene, sep = "\t", header = F, as.is = T)
colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
tg = tg[tg$type == 'cds',]
grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

fn = file.path(cfg$cdir, "../41_novseq/21.bed")
tn = read.table(fn, sep = "\t", header = F, as.is = T)
grn = with(tn, GRanges(seqnames = V1, ranges = IRanges(V2+1, end = V3)))

olens = intersect_basepair(grg, grn)
#ds = ddply(cbind(tg, olen = olens), .(fam), summarise, olen = sum(olen), alen = sum(end - beg + 1))
#ds = cbind(ds, prop = ds$olen / ds$alen, org = qname)
gb = group_by(cbind(tg, olen = olens), id)
dx = dplyr::summarise(gb, len = sum(end-beg+1), olen = sum(olen), pct = olen/len, fam = fam[1])

ds = ddply(dx, .(fam), summarise, cnt = length(pct), prop = sum(pct>=0.5)/cnt)


do = rbind(do, cbind(ds, org = qname))
}
to = ddply(do, .(fam), summarise, q25 = quantile(prop, 0.25), q50 = quantile(prop, 0.5), q75 = quantile(prop, 0.75))

fo = file.path(diro, "42.genefam.novseq.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)


fi = file.path(diro, "42.genefam.novseq.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)

ffam = file.path(diro, "41.gene.fams.tbl")
fams = read.table(ffam, header = F, sep = "\t", as.is = T)[,1]

tis = ti[ti$fam %in% fams,]
tis = tis[order(tis$q50, tis$q25, tis$q75, decreasing = T),]
tis$fam = factor(tis$fam, levels = tis$fam)
p3 = ggplot(tis) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', position = 'dodge', geom_params = list(width = 0.7, size = 0.3)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01)) +
  scale_y_continuous(name = 'Proportion novel (absent in HM101)', expand = c(0.002, 0.002)) +
  theme_bw() +
#  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = file.path(diro, "18_novseq_genefam.pdf")
ggsave(p3, filename = fp, width = 5, height = 4)


##### enrichment anlaysis (old: quantile of 15 accessions)
do = data.frame()
for (qname in qnames_all) {
#qname = qnames_all[14]
  cfg = cfgs[[qname]]

  tg = read.table(cfg$gene, sep = "\t", header = F, as.is = T)
  colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
  tg = tg[tg$type == 'cds',]
  grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

  tcnt = length(unique(tg$id))
  x1 = ddply(tg, .(fam), summarise, pct_tot_cnt = length(unique(id)) / tcnt)
  
  tbps = sum(tg$end - tg$beg + 1)
  y1 = ddply(tg, .(fam), summarise, pct_tot_bps = sum(end-beg+1)/tbps)
  
  fn = file.path(cfg$cdir, "../41_novseq/21.bed")
  tn = read.table(fn, sep = "\t", header = F, as.is = T)
  grn = with(tn, GRanges(seqnames = V1, ranges = IRanges(V2+1, end = V3)))

  olens = intersect_basepair(grg, grn)
  
  gb = group_by(cbind(tg, olen = olens), id)
  ds = dplyr::summarise(gb, len = sum(end-beg+1), olen = sum(olen), pct = olen/len, fam = fam[1])
  
  x2 = ddply(ds, .(fam), summarise, cnt_nov = sum(pct>=0.5))
  novcnt = sum(x2$cnt_nov)
  x2 = within(x2, {pct_nov_cnt = cnt_nov / novcnt})
  
  y2 = ddply(ds, .(fam), summarise, bps_nov = sum(olen))
  novbps = sum(y2$bps_nov)
  y2 = within(y2, {pct_nov_bps = bps_nov / novbps})
  
  cat(qname, tcnt, novcnt, "\n", sep = " ")

  ds1 = merge(x1, x2, by = 'fam')
  ds1 = within(ds1, {fold_cnt = pct_nov_cnt / pct_tot_cnt; org = qname})
  ds2 = merge(y1, y2, by = 'fam')
  ds2 = within(ds2, {fold_bps = pct_nov_bps / pct_tot_bps})
  ds = merge(ds1, ds2, by = 'fam')
  ds[order(ds$fold_cnt, decreasing = T), ][1:30,]
  
  do = rbind(do, ds)
}
to = ddply(do, .(fam), summarise, cnt_q25 = quantile(fold_cnt, 0.25), cnt_q50 = quantile(fold_cnt, 0.5), cnt_q75 = quantile(fold_cnt, 0.75), bps_q25 = quantile(fold_bps, 0.25), bps_q50 = quantile(fold_bps, 0.5), bps_q75 = quantile(fold_bps, 0.75))
to[order(to$cnt_q50, decreasing = T),][1:30,]

ffam = file.path(diro, "41.gene.fams.tbl")
fams = read.table(ffam, header = F, sep = "\t", as.is = T)[,1]

tos = to[to$fam %in% fams,]
tos = tos[order(tos$cnt_q50, tos$cnt_q25, tos$cnt_q75, decreasing = T),]
tos$fam = factor(tos$fam, levels = tos$fam)
tos = cbind(tos, enrich_cnt = 'y', enrich_bps = 'y', stringsAsFactors = F)
tos$enrich_cnt[tos$cnt_q50 < 1] = 'n'
tos$enrich_bps[tos$bps_q50 < 1] = 'n'
p1 = ggplot(tos) +
  geom_crossbar(aes(x = fam, y = cnt_q50, ymin = cnt_q25, ymax = cnt_q75, fill = enrich_cnt),
    stat = 'identity', position = 'dodge', geom_params = list(width = 0.7, size = 0.3)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01)) +
  scale_y_continuous(name = 'Fold change in novel gene pool (by gene count)', expand = c(0.02, 0.002)) +
  scale_fill_brewer(palette = "Set2", labels = c("Under-represented", "Over-represented")) +
  theme_bw() +  
  theme(legend.position = c(0.7,0.85), legend.direction = "vertical", legend.title = element_blank(), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.text = element_text(size = 8), legend.key.size = unit(0.8, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

p2 = ggplot(tos) +
  geom_crossbar(aes(x = fam, y = bps_q50, ymin = bps_q25, ymax = bps_q75, fill = enrich_bps),
    stat = 'identity', position = 'dodge', geom_params = list(width = 0.7, size = 0.3)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01)) +
  scale_y_continuous(name = 'Fold change in novel gene pool (by base pairs)', expand = c(0.02, 0.002)) +
  scale_fill_brewer(palette = "Set2", labels = c("Under-represented", "Over-represented")) +
  theme_bw() +  
  theme(legend.position = "none", legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 8), legend.key.size = unit(0.8, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
#  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))
  theme(axis.text.y = element_blank())

fp = file.path(diro, "18_nov_genefam.pdf")
numcol = 2
wds = c(4.5, 3.5)
pdf(file = fp, width = sum(wds), height = 5, bg = 'transparent')
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, numcol, width = wds)))

print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

dco = data.frame(x = rep(1,numcol), y = 1:numcol, lab = LETTERS[1:numcol])
for (i in 1:nrow(dco)) {
  x = dco$x[i]; y = dco$y[i]; lab = dco$lab[i]
  grid.text(lab, x = 0, y = unit(1, 'npc'), just = c('left', 'top'), gp = gpar(col = "black", fontface = 2, fontsize = 20),
    vp = viewport(layout.pos.row = x, layout.pos.col = y))
}
dev.off()

##### enrichment analysis (contribution to novel gene pool: cnt + seq)
fi = file.path(dirw, '32.global.tbl')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ddply(ti, .(org), summarise, len = sum(end - beg + 1))
ti = within(ti, {zchr = paste(org, chr, sep = "_")})
grn = with(ti, GRanges(seqnames = zchr, ranges = IRanges(beg, end = end)))
cat(sum(width(grn)), sum(width(reduce(grn))), "\n")

gb = group_by(ti, cid)
dcl = dplyr::summarise(gb, n_org = n(), 
#  orgs = paste(sort(unique(as.character(org))), collapse = "_"),
  size = sum(end - beg + 1))

# read in all genes
tg = data.frame()
for (qname in qnames_all) {
  cfg = cfgs[[qname]]
  tg1 = read.table(cfg$gene, sep = "\t", header = F, as.is = T)
  colnames(tg1) = c("chr", "beg", "end", "srd", "id", "type", "fam")
  tg1 = tg1[tg1$type == 'cds',]
  tg1 = within(tg1, {
    zchr = paste(qname, chr, sep = "_");
    zid = paste(qname, id, sep = "_")
  })
  tg = rbind(tg, tg1)
}
grg = with(tg, GRanges(seqnames = zchr, ranges = IRanges(beg, end = end)))
tlen = sum(tg$end - tg$beg + 1)
x1 = ddply(tg, .(fam), summarise, pct = sum(end-beg+1) / tlen)
gb = group_by(tg, zid)
x2 = dplyr::summarise(gb, fam = fam[1], glen = sum(end-beg+1))
x3 = ddply(x2, .(fam), summarise, pct = length(zid) / nrow(x2))

  t1 = ti[,c('zchr','beg','end','cid')]
  t2 = tg[,c('zchr','beg','end','zid','fam')]
  t1$beg = t1$beg - 1
  t2$beg = t2$beg - 1

  fbd1 = 'xtest1.bed'
  fbd2 = 'xtest2.bed'
  fres = 'xout.bed'
  options(scipen = 999)
  write.table(t1, fbd1, sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(t2, fbd2, sep = "\t", row.names = F, col.names = F, quote = F)
  options(scipen = 0)
  system(sprintf("intersectBed -wo -a %s -b %s > %s", fbd1, fbd2, fres))

  t3 = read.table(fres, sep = "\t", header = F, as.is = T)
  colnames(t3) = c('chr', 'beg1', 'end1', 'cid', 'chr2', 'beg2', 'end2', 'zid', 'fam', 'olen')
  system(sprintf("rm %s %s %s", fbd1, fbd2, fres))

## contribution to novel sequence content (bps)
to = within(t3, {
  beg = pmax(beg1, beg2);
  end = pmin(end1, end2);
  rm(beg1, end1, chr2, beg2, end2)
})
to = merge(to, dcl[,1:2], by = 'cid')
to = within(to, {polen = olen / n_org})

tlen = sum(to$polen)
y1 = ddply(to, .(fam), summarise, pct = sum(polen) / tlen)

tp = merge(x1, y1, by = 'fam')
colnames(tp) = c('fam', 'pct_total', 'pct_nov')
tp = within(tp, {fold = pct_nov / pct_total})
tp[order(tp$fold, decreasing = T), ][1:30,]

## contribution to novel gene pool by count (gene number)
to = within(t3, {
  beg = pmax(beg1, beg2);
  end = pmin(end1, end2);
  rm(beg1, end1, chr2, beg2, end2)
})
to = merge(to, dcl[,1:2], by = 'cid')
gb = group_by(to, cid, zid)
to = dplyr::summarise(gb, olen = sum(olen), n_org = n_org[1])
to = merge(to, x2, by = 'zid')
to = to[to$olen / to$glen >= 0.5,]

tcnt = nrow(to)#sum(1/to$n_org)
y2 = ddply(to, .(fam), summarise, pct = length(zid) / tcnt)#sum(1/n_org) / tcnt)

tp = merge(x3, y2, by = 'fam')
colnames(tp) = c('fam', 'pct_total', 'pct_nov')
tp = within(tp, {fold = pct_nov / pct_total})
tp[order(tp$fold, decreasing = T), ][1:30,]



ffam = file.path(diro, "41.gene.fams.tbl")
fams = read.table(ffam, header = F, sep = "\t", as.is = T)[,1]

tos = tp[tp$fam %in% fams,]
tos = tos[order(tos$fold, decreasing = T),]
tos$fam = factor(tos$fam, levels = tos$fam)
tos = cbind(tos, enrich = 'y', stringsAsFactors = F)
tos$enrich[tos$fold < 1] = 'n'
p3 = ggplot(tos) +
  geom_bar(aes(x = fam, y = fold, fill = enrich),
    stat = 'identity', geom_params = list(width = 0.7, size = 0.3)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01)) +
  scale_y_continuous(name = 'Fold change in novel gene pool (by sequence content)', expand = c(0.02, 0.002)) +
  scale_fill_brewer(palette = "Set2", labels = c("Under-represented", "Over-represented")) +
  theme_bw() +  
  theme(legend.position = "top", legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 8), legend.key.size = unit(0.8, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,0.5,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = file.path(diro, "18_nov_genefam_seq.pdf")
ggsave(p3, filename = fp, width = 5, height = 5)

