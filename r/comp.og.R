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
#write.table(tc, fo, sep = "\t", row.names = F, col.names = T, quote = F)

#### pan proteome revisited
fi = file.path(dirw, "05.clu/32.tbl")
ti = read.table(fi, sep = "\t", header = T, as.is = T)
x = strsplit(ti$id, "-")
ti = cbind(ti, org = sapply(x, "[", 1), gid = sapply(x, "[", 2))

ti = ti[ti$org %in% c(tname,qnames_12),]
ti = merge(ti, tc, by = c('org', 'gid'))

gb = group_by(ti, grp)
tr = dplyr::summarise(gb, size = length(unique(org)), org = org[1], fam = names(sort(table(cat2), decreasing = T))[1], rid = id[1])
ti = merge(ti, tr[,c('grp','size')], by = 'grp')

# output unknown cluster/ids for blastnr
ids_unk = tr$rid[tr$fam == 'Unknown']
fo = file.path(dirw, "31.unk.txt")
#write(ids_unk, fo, sep = "\n")
#seqret.pl -d 02.fas -b 31.unk.txt -o 32.unk.fas

tr$fam[tr$fam %in% c("NB-ARC", "TIR", 'CC-NBS-LRR','TIR-NBS-LRR')] = "NBS-LRR"
fams = c("NBS-LRR", "LRR", "F-box", "LRR-RLK", "NCR", "TE", "Unknown")
tr$fam[! tr$fam %in% fams] = 'Pfam-Miscellaneous'
tr$fam = factor(tr$fam, levels = c(fams, 'Pfam-Miscellaneous'))

table(tr$size)
table(tr$org[tr$size==1])

grps = tr$grp[tr$size == 1]
ts = tr[tr$grp %in% grps,]

y = table(ts$fam)
y[order(y, decreasing = T)][1:20]/sum(y)

brks = as.character(1:13)

famss = c('Pfam-Miscellaneous', 'F-box', 'LRR-RLK', 'NCR', 'NBS-LRR', 'Unknown')
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
  scale_y_continuous(name = '# Clusters') +
  facet_wrap(~ fam, scales = 'free', nrow = 2) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 90, hjust = 1))

fp = sprintf("%s/41.afs.pdf", dirw)
ggsave(p1, filename = fp, width = 7, height = 5)


##### double check unknown proteins
fq = file.path(Sys.getenv("genome"), "HM101", "augustus", "31.gtb")
tq = read.table(fq, sep = "\t", header = T, as.is = T)

idxs = which(grepl("\\[\\w+\\]", tq$note))
m = regexpr("\\[(\\w+)\\]", tq$note, perl = TRUE)
x = regmatches(tq$note, m)
stopifnot(length(x) == length(idxs))
tqs = data.frame(id = tq$id[idxs], qual = x, note = tq$note[idxs])
ids_hc = paste(tname, tqs$id[tqs$qual == '[HC]'], sep = "-")
ids_lc = paste(tname, tqs$id[tqs$qual == '[LC]'], sep = "-")

x = strsplit(tr$rid, "-")
gids = sapply(x, "[", 2)
#table(tqs$qual[tqs$id %in% gids[tr$org == 'HM101' & tr$size >1 & tr$fam != 'Unknown']])


ds = data.frame()
for (qname in qnames_12) {
  fq = file.path(Sys.getenv("genome"), qname, "augustus", "31.gtb")
  tq = read.table(fq, sep = "\t", header = T, as.is = T)
  dss = tq[,c('id','note')]
  colnames(dss)[2] = 'qual'
  dss = cbind(dss, rid = paste(qname, dss$id, sep = "-"), stringsAsFactors = F)
  ds = rbind(ds, dss)
}
#hist(ds$qual[ds$rid %in% ti$id[ti$size == 1 & ti$cat2 == 'Unknown']])
#hist(ds$qual[ds$rid %in% ti$id[ti$size == 1 & ti$cat2 != 'Unknown']])
#hist(ds$qual[ds$rid %in% ti$id[ti$size == 13 & ti$cat2 == 'Unknown']])
ids_hc = c(ids_hc, ds$rid[ds$qual >= 0.9])
ids_lc = c(ids_lc, ds$rid[ds$qual < 0.9])


do = data.frame()
for (qname in c("HM056")) {
  fx = file.path(Sys.getenv("misc2"), "rnaseq/mt/31_cufflinks", qname, "isoforms.fpkm_tracking")
  tx = read.table(fx, sep = "\t", header = T, as.is = T)
  dos = data.frame(rid = paste(qname, tx$tracking_id, sep = "-"), fpkm = tx$FPKM, stringsAsFactors = F)
  do = rbind(do, dos)
}
do$fpkm[do$fpkm > 0] = 1
hist(do$fpkm[do$rid %in% ti$id[ti$size == 1 & ti$cat2 != 'Unknown']])
hist(do$fpkm[do$rid %in% ti$id[ti$size == 1 & ti$cat2 == 'Unknown']])
hist(do$fpkm[do$rid %in% ti$id[ti$size > 12 & ti$cat2 == 'Unknown']])


t.unk = tr[tr$fam == 'Pfam-Miscellaneous',]
t.unk = tr[tr$fam == 'Unknown',]
t.unk.hc = t.unk[t.unk$rid %in% ids_hc,]
t.unk.lc = t.unk[t.unk$rid %in% ids_lc,]

brks = as.character(1:13)

x = table(t.unk$size)
do1 = data.frame(lab = 'All Unknown Proteins', size = brks, cnt = as.numeric(x[brks]), stringsAsFactors = F)
x = table(t.unk.hc$size)
do2 = data.frame(lab = 'High Confidence', size = brks, cnt = as.numeric(x[brks]), stringsAsFactors = F)
x = table(t.unk.lc$size)
do3 = data.frame(lab = 'Low Confidence', size = brks, cnt = as.numeric(x[brks]), stringsAsFactors = F)
do = rbind(do2, do3)
do$size = factor(do$size, levels = brks)

p1 = ggplot(do, aes(x = size, y = cnt, fill = lab)) +
  geom_bar(stat = 'identity', position = 'stack', geom_params=list(width = 0.8)) +
  scale_x_discrete(name = '', breaks = brks) +
  scale_y_continuous(name = '# Clusters') +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = c(0.6, 0.8), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.6, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_blank(), legend.text = element_text(size = 9, angle = 0)) +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 90, hjust = 1))

fp = sprintf("%s/45.unk.pdf", dirw)
ggsave(p1, filename = fp, width = 4, height = 4)


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


tr$fam[tr$fam %in% c("NB-ARC", "TIR", 'CC-NBS-LRR','TIR-NBS-LRR')] = "NBS-LRR"
fams = c("NCR", "NBS-LRR", "LRR", "F-box", "LRR-RLK", "TE", "Unknown")
tr$fam[! tr$fam %in% fams] = 'Pfam-Miscellaneous'
tr$fam = factor(tr$fam, levels = c(fams, 'Pfam-Miscellaneous'))

x = table(tr$fam)
tp = data.frame(fam = names(x), num = as.numeric(x), stringsAsFactors = T)
tp = cbind(tp, prop = tp$num / sum(tp$num))
labs = sprintf("%s (%d : %.01f%%)", tp$fam, tp$num, tp$prop * 100)
tp = cbind(tp, label = labs)

fo = file.path(dirw, "51.novel.fam.tbl")
write.table(tp, fo, sep = "\t", row.names = F, col.names = T, quote = F)

fp = file.path(dirw, '51.novel.fam.pdf')
pdf(file = fp, width = 8, height = 6)
pie(tp$num, labels = labs, lwd = 1, cex = 0.7, col = brewer.pal(8, 'Set3'), main = sprintf("Composition of accession-specific gene pool (%d)", nrow(tr)))
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