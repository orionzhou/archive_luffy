require(GenomicRanges)
require(plyr)
require(rtracklayer)
require(Rsamtools)
require(ggplot2)
require(grid)

orgs = c(
  "HM058", "HM125", "HM056", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340"
)
chrs = sprintf("chr%s", 1:8)


##### Mapping-based approach

dirr = file.path(Sys.getenv('genome'), "HM101")
fa = file.path(dirr, '15.sizes')
fg = file.path(dirr, '16.gap.bed')
ta = read.table(fa, sep = '\t', header = F, as.is = T)
tg = read.table(fg, sep = '\t', header = F, as.is = T)
ta = ta[ta$V1 %in% chrs,]
tg = tg[tg$V1 %in% chrs,]
ga = GRanges(seqnames = ta$V1, ranges = IRanges(1, end = ta$V2))
gg = GRanges(seqnames = tg$V1, ranges = IRanges(tg$V2+1, end = tg$V3))
gr = setdiff(ga, gg)

ff = file.path(dirr, "51.tbl")
tf = read.table(ff, sep = '\t', header = F, as.is = T)
colnames(tf) = c("chr", "beg", "end", "srd", "id", "type", "cat")

fc = file.path(dirr, "51.merged.tbl")
tc = read.table(fc, sep = '\t', header = T, as.is = T)
gc = GRanges(seqnames = tc$chr, ranges = IRanges(tc$beg, end = tc$end))

dirm = file.path(Sys.getenv("misc3"), "hapmap", "12_ncgr")

### generate coverage BED file (takes long time - run with caution)
for (org in orgs) {
fcov = sprintf("%s/35_cov/%s.bw", dirm, org)

bw = import(fcov, which = gr, asRangedData = F)
chrs = seqnames(bw)
poss = start(bw)
idxs = bw$score >= 1
rm(bw)
gm = GRanges(seqnames = chrs[idxs], ranges = IRanges(poss[idxs], 
  end = poss[idxs]))
gm = reduce(gm)
rm(chrs, poss, idxs)

tm = data.frame(chr = seqnames(gm), beg = start(gm) - 1, end = end(gm))
fo = sprintf("%s/38_covered/%s.bed", dirm, org)
write.table(tm, fo, row.names = F, col.names = F, sep = "\t", quote = F)
}

### SNP density
tr = data.frame()

for (org in orgs) {
#org = "HM004"

fs = file.path(dirm, "42_vnt", org, "snp")
ts = read.table(fs, sep = "\t", header = F, as.is = T)
colnames(ts) = c("chr", "pos", "ref", "alt")
gs = GRanges(seqnames = ts$chr, ranges = IRanges(ts$pos, end = ts$pos))

fv = sprintf("%s/38_covered/%s.bed", dirm, org)
tv = read.table(fv, sep = "\t", header = F, as.is = T)
gv = GRanges(seqnames = tv$V1, ranges = IRanges(tv$V2+1, end = tv$V3))

for (type in unique(tc$type)) {
  gcs = gc[tc$type == type]
  gsi = intersect(gs, gcs)
  nsnp = sum(width(gsi))
  
  gvi = intersect(gv, gcs)
  bcov = sum(width(gvi))
  trs = data.frame(org = org, type = type, bcov = bcov, nsnp = nsnp)
  tr = rbind(tr, trs)
}
}

tr = cbind(tr, subrate = tr$nsnp / tr$bcov)

##### Assembly-based approach
### SNP density
ta = data.frame()

for (org in orgs) {
#org = "HM034"
dira = sprintf("%s/%s_HM101/23_blat/31.9", Sys.getenv('misc3'), org)

fs = file.path(dira, "snp")
ts = read.table(fs, sep = "\t", header = F, as.is = T)
colnames(ts) = c("chr", "pos", "ref", "alt", "qid", "qpos", "cid", "lev")
ts = ts[ts$lev <= 2,]
gs = GRanges(seqnames = ts$chr, ranges = IRanges(ts$pos, end = ts$pos))

fv = file.path(dira, "gax")
tv = read.table(fv, sep = "\t", header = F, as.is = T)
colnames(tv) = c("chr", "beg", "end", "srd", "qid", "qbeg", "qend", "qsrd",
  "cid", "lev")
tv = tv[tv$lev <= 2,]
gv = GRanges(seqnames = tv$chr, ranges = IRanges(tv$beg, end = tv$end))

for (type in unique(tc$type)) {
  gcs = gc[tc$type == type]
  gsi = intersect(gs, gcs)
  nsnp = sum(width(gsi))
  
  gvi = intersect(gv, gcs)
  bcov = sum(width(gvi))
  tas = data.frame(org = org, type = type, bcov = bcov, nsnp = nsnp)
  ta = rbind(ta, tas)
}
}

ta = cbind(ta, subrate = ta$nsnp / ta$bcov)


##### compare SNP density estimates from two sources
to = rbind(cbind(method = 'mapping', tr), cbind(method = 'assembly', ta))
fo = file.path(Sys.getenv("misc3"), "compstat", "snp.pdf")

to$type = factor(to$type, levels = c('CDS', 'Intron', 'UTR', 'Intergenic'))
to$method = factor(to$method, levels = c("mapping", "assembly"))
p = ggplot(to) +
  geom_bar(mapping = aes(x = type, y = subrate, fill = method), 
    stat = 'identity', position = 'dodge', geom_params=list(width = 0.5)) +
  scale_fill_brewer(palette='Set2', name = '',
#    breaks = c('assembly', 'mapping'), 
    labels = c("Mapping-based", "Assembly-based")) +
  scale_x_discrete(name = '') +
  scale_y_continuous(name = 'SNP density (/bp)') +
  facet_wrap( ~ org, nrow = 5) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(0,1,1,0), "lines")) +
  theme(axis.text.x = element_text(size = 8, colour = "grey", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0)) +
  theme(panel.margin = unit(0.6, 'lines'))
ggsave(p, filename = fo, width = 6, height = 6)


##### compare SNPs/Indels called by two approaches
do = data.frame()
for (org in orgs) {
#org = "HM034"
  dira = sprintf("%s/%s_HM101/23_blat/31.9", Sys.getenv('misc3'), org)
  fa = file.path(dira, "vnt.tbl")
  ta = read.table(fa, sep = "\t", header = F, as.is = T)
  colnames(ta) = c("chr", "pos", "ref", "alt", "score")

  dirm = file.path(Sys.getenv("misc3"), "hapmap", "12_ncgr", "44_tbl")
  fm = sprintf("%s/%s.tbl", dirm, org)
  tm = read.table(fm, sep = "\t", header = F, as.is = T)
  colnames(tm) = c("chr", "pos", "ref", "alt", "score")

  tms = tm[tm$chr %in% chrs & nchar(tm$ref) == 1 & nchar(tm$alt) == 1, ]
  tas = ta[ta$chr %in% chrs & nchar(ta$ref) == 1 & nchar(ta$alt) == 1, ]
  tc = merge(tas[,1:4], tms, by = c('chr', 'pos'), all = T)
  
  t_ovl = tc[!is.na(tc$alt.x) & !is.na(tc$alt.y),]
  t_anm = tc[!is.na(tc$alt.x) & is.na(tc$alt.y),]
  t_mna = tc[is.na(tc$alt.x) & !is.na(tc$alt.y),]
  n_con = sum(t_ovl$alt.x == t_ovl$alt.y)
  n_dis = sum(t_ovl$alt.x != t_ovl$alt.y)

  d1 = data.frame(org = org, vnt_type = 'snp', type = c("Assembly-Only", 
    "Assembly+Mapping:Concordant", "Assembly+Mapping:Discordant", 
    "Mapping-Only"), cnt = c(nrow(t_anm), n_con, n_dis, nrow(t_mna)))
  
  tms = tm[tm$chr %in% chrs & (nchar(tm$ref) != 1 | nchar(tm$alt) != 1), ]
  tas = ta[ta$chr %in% chrs & (nchar(ta$ref) != 1 | nchar(ta$alt) != 1), ]
  tc = merge(tas[,1:4], tms, by = c('chr', 'pos'), all = T)
  
  t_ovl = tc[!is.na(tc$alt.x) & !is.na(tc$alt.y),]
  t_anm = tc[!is.na(tc$alt.x) & is.na(tc$alt.y),]
  t_mna = tc[is.na(tc$alt.x) & !is.na(tc$alt.y),]
  n_con = sum(t_ovl$ref.x == t_ovl$ref.y & t_ovl$alt.x == t_ovl$alt.y)
  n_dis = sum(t_ovl$ref.x != t_ovl$ref.y | t_ovl$alt.x != t_ovl$alt.y)

  d2 = data.frame(org = org, vnt_type = 'indel', type = c("Assembly-Only", 
    "Assembly+Mapping:Concordant", "Assembly+Mapping:Discordant", 
    "Mapping-Only"), cnt = c(nrow(t_anm), n_con, n_dis, nrow(t_mna)))
  do = rbind(do, d1, d2)
}

fo = file.path(Sys.getenv("misc3"), "compstat", "comp_vnt.pdf")

do$org = factor(do$org, levels = orgs)
#do$type = factor(do$type, levels = c('CDS', 'Intron', 'UTR', 'Intergenic'))
p = ggplot(do) +
  geom_bar(mapping = aes(x = org, y = cnt, fill = type), 
    stat = 'identity', position = 'stack', geom_params=list(width = 0.5)) +
  scale_fill_manual(values = c('skyblue1', 'firebrick1', 'orchid', 'palegreen'), name = '', guide = guide_legend(nrow = 2, byrow = T, label.position = "right", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  scale_x_discrete(name = '') +
  scale_y_continuous(name = 'Variant #') +
  facet_wrap( ~ vnt_type, ncol = 1) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(0,1,1,0), "lines")) +
  theme(axis.text.x = element_text(size = 8, colour = "grey", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0))
ggsave(p, filename = fo, width = 6, height = 6)


##### show variant quality score distribution of mapping/assembly
org = "HM034"
  dira = sprintf("%s/%s_HM101/23_blat/31.9", Sys.getenv('misc3'), org)
  fa = file.path(dira, "vnt.tbl")
  ta = read.table(fa, sep = "\t", header = F, as.is = T)
  colnames(ta) = c("chr", "pos", "ref", "alt", "score")

  dirm = file.path(Sys.getenv("misc3"), "hapmap", "12_ncgr", "44_tbl")
  fm = sprintf("%s/%s.tbl", dirm, org)
  tm = read.table(fm, sep = "\t", header = F, as.is = T)
  colnames(tm) = c("chr", "pos", "ref", "alt", "score")

  tms = tm[tm$chr %in% chrs & nchar(tm$ref) == 1 & nchar(tm$alt) == 1, ]
  tas = ta[ta$chr %in% chrs & nchar(ta$ref) == 1 & nchar(ta$alt) == 1, ]
  tc = merge(tas[,1:4], tms, by = c('chr', 'pos'), all = T)
  
  t_ovl = tc[!is.na(tc$alt.x) & !is.na(tc$alt.y),]
  t_anm = tc[!is.na(tc$alt.x) & is.na(tc$alt.y),]
  t_mna = tc[is.na(tc$alt.x) & !is.na(tc$alt.y),]

to = rbind(cbind(t_ovl, type = 'Assembly+Mapping'), 
  cbind(t_mna, type = 'Mapping-Only'))

fo = sprintf("%s/compstat/comp_snp_%s.pdf", Sys.getenv("misc3"), org)
p = ggplot(to) +
  geom_histogram(mapping = aes(x = score, fill = type), 
    position = "dodge", geom_params=list(width = 0.5, alpha = 0.8)) +
  scale_fill_manual(values = c('skyblue1', 'firebrick1', 'orchid', 'palegreen'), name = '', guide = guide_legend(nrow = 1, byrow = T, label.position = "right", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  scale_x_continuous(name = 'variant score', limits = c(0, 60000)) +
  scale_y_continuous(name = 'count') +
  theme_bw() +
  theme(legend.position = "top", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(0,1,1,0), "lines")) +
  theme(axis.text.x = element_text(size = 8, colour = "grey", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0))
ggsave(p, filename = fo, width = 6, height = 6)


##### compare InDels
chrs = sprintf("chr%s", 1:8)

do = data.frame()
for (org in c("HM058", "HM060", "HM095", "HM034")) {
#org = "HM034"
  dira = sprintf("%s/%s_HM101/23_blat/31.9", Sys.getenv('misc3'), org)
  fa = file.path(dira, "vnt.tbl")
  ta = read.table(fa, sep = "\t", header = F, as.is = T)
  colnames(ta) = c("chr", "pos", "ref", "alt", "score")

  dirm = file.path(Sys.getenv("misc3"), "hapmap", "12_ncgr", "44_tbl")
  fm = sprintf("%s/%s.tbl", dirm, org)
  tm = read.table(fm, sep = "\t", header = F, as.is = T)
  colnames(tm) = c("chr", "pos", "ref", "alt", "score")

  tms = tm[tm$chr %in% chrs & (nchar(tm$ref) != 1 | nchar(tm$alt) != 1), ]
  tas = ta[ta$chr %in% chrs & (nchar(ta$ref) != 1 | nchar(ta$alt) != 1), ]
#  tc = merge(tas[,1:4], tms, by = c('chr', 'pos'), all = T)
  
  x = c(-49:-1, 1:49)
  
  tb = table(nchar(tms$alt) - nchar(tms$ref))
  y = tb[as.character(x)]
  y[is.na(y)] = 0
  df1 = data.frame(type = 'Mapping-based', len = x, cnt = y)

  tb = table(nchar(tas$alt) - nchar(tas$ref))
  y = tb[as.character(x)]
  y[is.na(y)] = 0
  df2 = data.frame(type = 'Assembly-based', len = x, cnt = y)

  dos = cbind(rbind(df1, df2), org = org)
  do = rbind(do, dos)
}
#  do$cnt = log(do$cnt)

fo = sprintf("%s/compstat/comp_idm_size.pdf", Sys.getenv("misc3"))
p = ggplot(do) +
  geom_bar(mapping = aes(x = len, y = cnt, fill = type), 
    stat = 'identity', position = 'dodge', 
    geom_params=list(width = 0.8, alpha = 0.8)) +
  scale_fill_brewer(palette='Set1', name = '', guide = guide_legend(nrow = 1, byrow = T, label.position = "right", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  scale_x_continuous(name = 'Indel Size (bp)') +
  scale_y_continuous(name = '# events (log)') +
  facet_wrap( ~ org, nrow = 2) +
  theme_bw() +
  theme(legend.position = "top", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(0,1,1,0), "lines")) +
  theme(axis.text.x = element_text(size = 8, colour = "grey", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0))
ggsave(p, filename = fo, width = 6, height = 6)

