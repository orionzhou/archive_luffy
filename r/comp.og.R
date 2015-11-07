require(plyr)
require(dplyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(reshape2)
require(RColorBrewer)
require(ape)
require(gridBase)
require(colorRamps)
source("Location.R")
source("comp.fun.R")

dirw = file.path(Sys.getenv("misc3"), "comp.og")

##### dump all proteins
tc = data.frame()
for (qname in c(tname, qnames_15)) {
  gdir = cfgs[[qname]]$gdir

  fg = file.path(gdir, "51.gtb")
  tg = read.table(fg, sep = "\t", header = T, as.is = T)

  tc = rbind(tc, data.frame(org = qname, gid = tg$id, cat2 = tg$cat2, cat3 = tg$cat3, stringsAsFactors = F))
}
head(tc)
table(tc$org)

fo = file.path(dirw, "01.gid.tbl")
write.table(tc, fo, sep = "\t", row.names = F, col.names = T, quote = F)

#### pan proteome revisited
fi = file.path(dirw, "05.clu/32.tbl")
ti = read.table(fi, sep = "\t", header = T, as.is = T)
x = strsplit(ti$id, "-")
ti = cbind(ti, org = sapply(x, "[", 1), gid = sapply(x, "[", 2))

ti = ti[ti$org %in% c(tname,qnames_12),]
ti = merge(ti, tc, by = c('org', 'gid'))

gb = group_by(ti, grp)
tr = dplyr::summarise(gb, size = length(unique(org)), org = org[1], fam = names(sort(table(cat2), decreasing = T))[1])

tr$fam[tr$fam %in% c('CC-NBS-LRR','TIR-NBS-LRR')] = "NBS-LRR"
tr$fam[tr$fam %in% c("NB-ARC", "TIR")] = "NB-ARC / TIR"
fams = c("NBS-LRR", "NBS-LRR", "NB-ARC / TIR", "LRR", "F-box", "LRR-RLK", "NCR", "TE", "Unknown")
tr$fam[! tr$fam %in% fams] = 'Pfam-Miscellaneous'
tr$fam = factor(tr$fam, levels = c(fams, 'Pfam-Miscellaneous'))

table(tr$size)
table(tr$org[tr$size==1])

grps = tr$grp[tr$size == 1]
ts = tr[tr$grp %in% grps,]

y = table(ts$fam)
y[order(y, decreasing = T)][1:20]/sum(y)

brks = as.character(1:13)

famss = c('Pfam-Miscellaneous', 'NCR', 'NBS-LRR', 'Unknown')
do = data.frame()
for (fam in famss) {
  x = table(tr$size[tr$fam == fam])
  dos = data.frame(fam = fam, size = brks, cnt = as.numeric(x[brks]), stringsAsFactors = F)
  do = rbind(do, dos)
}
do$size = factor(do$size, levels = brks)
do$fam = factor(do$fam, levels = famss)

p1 = ggplot(do, aes(x = size, y = cnt)) +
  geom_bar(stat = 'identity', geom_params=list(width = 0.8)) +
  scale_x_discrete(name = '', breaks = brks) +
  scale_y_continuous(name = '# Groups') +
  facet_wrap(~ fam, scales = 'free', nrow = 2) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 90, hjust = 1))

fp = sprintf("%s/41.afs.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 5)



tab1 = table(tr$size)
tab1 = tab1[names(tab1) != 1]
dt1 = data.frame(norg = as.numeric(names(tab1)), cnt = as.numeric(tab1), org = 'mixed', stringsAsFactors = F)

tab2 = table(tr$org[tr$size == 1])
tab2 = tab2[tab2 > 0]
dt2 = data.frame(norg = 1, cnt = as.numeric(tab2), org = names(tab2), stringsAsFactors = F)

to = rbind(dt1, dt2)

cols = c(brewer.pal(12, 'Set3'), brewer.pal(3, 'Set1')[1], 'gray30')
labs = orgs

to$org = factor(to$org, levels = c(orgs, 'mixed'))
to$norg = factor(to$norg, levels = sort(as.numeric(unique(to$norg))))
p1 = ggplot(to, aes(x = norg, y = cnt, fill = org, order = plyr:::desc(org))) +
  geom_bar(stat = 'identity', position = "stack", geom_params=list(width = 0.5)) +
  scale_fill_manual(name = "Accession-Specific", breaks = labs, labels = labs, values = cols, guide = guide_legend(ncol = 1, byrow = F, label.position = "right", direction = "vertical", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  scale_x_discrete(name = '# Sharing Accession') +
  scale_y_continuous(name = '# Gene Models', expand = c(0, 0), limits = c(0, 30100)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(legend.position = c(0.3, 0.7), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.6, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_text(size = 8, angle = 0), legend.text = element_text(size = 8, angle = 0)) +
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 90, hjust = 1))

fp = sprintf("%s/40.pan.proteome.afs.pdf", dirw)
ggsave(p1, filename = fp, width = 5, height = 5)


##### identify all novel / accession-specific genes
fi = file.path(Sys.getenv('misc3'), 'comp.panseq', '32.global.tbl')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

gb = group_by(ti, cid)
dcl = dplyr::summarise(gb, n_org = n(), size = sum(end - beg + 1))
cids = dcl$cid[dcl$n_org == 1]

do = data.frame()
for (qname in qnames_15) {
  cfg = cfgs[[qname]]

  tg = read.table(cfg$gene, sep = "\t", header = F, as.is = T)
  colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
  tg = tg[tg$type == 'cds',]
  grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

  fn = file.path(cfg$cdir, "../41_novseq/21.bed")
  tn = read.table(fn, sep = "\t", header = F, as.is = T)
  grn = with(tn, GRanges(seqnames = V1, ranges = IRanges(V2+1, end = V3)))
#  tn = ti[ti$cid %in% cids & ti$org == qname,]
#  grn = with(tn, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

  olens = intersect_basepair(grg, grn)
  #ds = ddply(cbind(tg, olen = olens), .(fam), summarise, olen = sum(olen), alen = sum(end - beg + 1))
  #ds = cbind(ds, prop = ds$olen / ds$alen, org = qname)
  gb = group_by(cbind(tg, olen = olens), id)
  dx = dplyr::summarise(gb, len = sum(end-beg+1), olen = sum(olen), pct = olen/len, fam = fam[1])

  dy = dx[dx$pct >= 0.5,]
  do = rbind(do, cbind(dy, org = qname))
}
ids_nov = sprintf("%s-%s", do$org, do$id)

x = ddply(do, .(fam), summarise, cnt = length(fam), prop = cnt / nrow(do))
x[order(x$prop, decreasing = T), ][1:30,]

fn = file.path(Sys.getenv("misc3"), "comp.ortho", "06.no.ortho.tbl")
tn = read.table(fn, sep = "\t", header = T, as.is = T)
ids2 = sprintf("%s-%s", tn$org, tn$gid)

sum(ids_nov %in% ids2)

to = data.frame(org = do$org, gid = do$id, stringsAsFactors = F)
fo = file.path(dirw, "11.novel.gid.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)


##### run seq.cluster.py and characterize novel gene pool composition
fi = file.path(dirw, "15.clu/32.tbl")
ti = read.table(fi, sep = "\t", header = T, as.is = T)
x = strsplit(ti$id, "-")
ti = cbind(ti, org = sapply(x, "[", 1), gid = sapply(x, "[", 2))

fd = file.path(dirw, "11.novel.gid.tbl")
td = read.table(fd, sep = "\t", header = T, as.is = T)

tt = merge(td, ti, by = c("org", "gid"), all.x = T)
maxgrp = max(tt$grp, na.rm = T)
nna = sum(is.na(tt$grp))
tt$grp[is.na(tt$grp)] = seq(maxgrp+1, maxgrp+nna)

gb = group_by(tt, grp)
x = dplyr::summarise(gb, size = length(unique(org)))
table(x$size)
grps = x$grp[x$size == 1]

ti2 = merge(tt[tt$grp %in% grps,], tc, by = c('org', 'gid'))
gb = group_by(ti2, grp)
tr = summarise(gb, fam = names(sort(table(cat2), decreasing = T))[1], sfam = cat3[which(cat2==fam)[1]])

tr$fam[tr$fam %in% c("NB-ARC", "TIR")] = "NB-ARC / TIR"
fams = c("CC-NBS-LRR", "TIR-NBS-LRR", "NB-ARC / TIR", "LRR", "F-box", "LRR-RLK", "NCR", "TE", "Unknown")
tr$fam[! tr$fam %in% fams] = 'Pfam-Miscellaneous'
tr$fam = factor(tr$fam, levels = c(fams, 'Pfam-Miscellaneous'))

x = table(tr$fam)
tp = data.frame(fam = names(x), num = as.numeric(x), stringsAsFactors = T)
tp = cbind(tp, prop = tp$num / sum(tp$num))
labs = sprintf("%s (%d : %.01f%%)", tp$fam, tp$num, tp$prop * 100)

fo = file.path(dirw, '51.novel.fam.pdf')
pdf(file = fo, width = 8, height = 6)
pie(tp$num, labels = labs, lwd = 1, cex = 0.7, main = sprintf("Composition of accession-specific gene pool (%d)", nrow(tr)))
dev.off()


fd = file.path(dirw, "01.gid.tbl.bak")
td = read.table(fd, sep = "\t", header = T, as.is = T)
ids_old = td$gid[td$org == 'HM056']
ids_new = tc$gid[tc$org == "HM056"]
ids_new[!ids_new %in% ids_old]

t1 = read.table(file.path(dirw, "11.novel.gid.tbl"), sep = "\t", header = T, as.is = T)
t2 = read.table(file.path(dirw, "16.len.tbl"), sep = "\t", header = F, as.is = T)
t3 = read.table(file.path(dirw, "16.tbl"), sep = "\t", header = F, as.is = T)

ids_all = paste(t1$org, t1$gid, sep = "-")
ids_cen = t2$V1
length(ids_all)
length(ids_cen)
sum(ids_cen %in% ids_all)

ids_hit = t3$V1
length(ids_hit)
length(unique(ids_hit))

length(unique(t3$V2))
sum(unique(t3$V2) %in% ids_cen)

ids_lef = ids_all[!ids_all %in% c(ids_hit, ids_cen)]