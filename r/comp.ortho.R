require(plyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(RColorBrewer)
require(ape)
require(gridBase)
require(colorRamps)

diro = file.path(Sys.getenv("misc3"), "comp.ortho")

tname = "HM101"
qnames = c(
  "HM058", "HM125", "HM056", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340"
)

##### create raw ortholog groups for 16 accessions & write un-ortholog seqs
f_tgene = file.path(Sys.getenv("genome"), tname, "51.gtb")
tgene = read.table(f_tgene, header = T, sep = "\t", as.is = T)[,c(1,16:17)]

to = tgene[,1:2]
colnames(to) = c(tname, 'cat')
i = 1
for (qname in qnames) {
  fi = sprintf("%s/%s_%s/51_ortho/05.score.tbl", Sys.getenv("misc3"), qname, tname)
  ti = read.table(fi, sep = "\t", as.is = T, header = T)
  ti = cbind(ti, tcov = (ti$mat+ti$mis)/ti$tlen, qcov = (ti$mat+ti$mis)/ti$qlen)
  idxs = (ti$qcov>=0.5 | ti$tcov>=0.5)
  n_pass = sum(idxs)
  n_fail = nrow(ti) - n_pass
  cat(sprintf("%s: %5d filtered, %5d passed\n", qname, n_fail, n_pass))
  tis = ti[idxs,c('tid','qid')]
  
  to = merge(to, tis, by.x = tname, by.y = 'tid', all.x = T)
  colnames(to)[2+i] = qname
  i = i + 1
}

n_org = apply(to, 1, function(z) sum(!is.na(z[c(-1,-2)])))
tt = cbind(to, n_org = n_org)
ft = file.path(diro, "01.ortho.tbl")
write.table(tt, ft, sep = "\t", row.names = F, col.names = T, quote = F)


ids = tt[tt$n_org != length(qnames), tname]
tp = data.frame(org = tname, id = ids)
cat(sprintf("%s: %5d added\n", tname, length(ids)))

for (qname in qnames) {
  f_qgene = file.path(Sys.getenv("genome"), qname, "51.gtb")
  qgenes = read.table(f_qgene, header = T, sep = "\t", as.is = T)[1]$id

  ids_o = tt[!is.na(tt[qname]), qname]
  ids = qgenes[!qgenes %in% ids_o]
  
  tps = data.frame(org = qname, id = ids)
  tp = rbind(tp, tps)
  cat(sprintf("%s: %5d added\n", qname, length(ids)))
}
fp = file.path(diro, "06.no.ortho.tbl")
write.table(tp, fp, sep = "\t", row.names = F, col.names = T, quote = F)


##### generate final score matrix of ortho groups
fi = file.path(diro, "31.ortho.tbl")
ti = read.table(fi, sep = "\t", header = T, stringsAsFactors = F)
ti[is.na(ti)] = ''

ma = matrix(NA, nrow(ti), ncol(ti))
colnames(ma) = colnames(ti)
ma[ti[,tname] != '', tname] = 1
ma[ti[,tname] == '', tname] = 0

for (qname in qnames) { 
  fr = sprintf("%s/%s_%s/51_ortho/05.score.tbl", Sys.getenv("misc3"), qname, tname)
  tr = read.table(fr, header = T, sep = "\t", stringsAsFactors = F)
  tr = cbind(tr, score = tr$mat / tr$len)
  
  tis = cbind(idx = 1:nrow(ti), ti[,c(tname, qname)], stringsAsFactors = F)
  tx = merge(tis, tr[,c('tid','qid','score')], 
    by.x = c(tname, qname), by.y = c('tid', 'qid'), all.x = T)
  ma[,qname] = tx$score[order(tx$idx)]
  ma[ti[,qname] == '', qname] = 0
  ma[ti[,tname] == '' & ti[,qname] != '', qname] = 1
}
sum(is.na(ma))

fo = file.path(diro, "35.ortho.score.tbl")
write.table(ma, fo, sep = "\t", row.names = F, col.names = T, quote = F)


##### CRP ortho-map HC
fi = file.path(diro, "33.ortho.cat.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti[is.na(ti)] = ''

fr = file.path(diro, "35.ortho.score.tbl")
tr = read.table(fr, sep = "\t", header = T, stringsAsFactors = F)

idxs = ti$cat2 == "CRP"
tis = ti[idxs,]
trs = tr[idxs,]
idxs = order(tis$cat2, tis$cat3, tis[,tname])
ti = tis[idxs,]
tr = trs[idxs,]
ti = cbind(ti, x = 1:nrow(ti), stringsAsFactors = F)

### hclust
cordist <- function(x) as.dist(1-cor(t(x), method="pearson"))
r.dist <- cordist(t(tr))
r.hc <- hclust(r.dist, method='ward.D')
r.dendro <- as.dendrogram(r.hc)


ff = '/home/youngn/zhoup/Data/misc2/genefam/crp.info'
tf = read.table(ff, sep = "\t", as.is = T, header = T)
tf = merge(ti[,c('cat3','x')], tf, by.x = 'cat3', by.y = 'id')
tx = ddply(tf, .(cat), summarise, bed = min(x), end = max(x), mid = as.integer(mean(x)))

tp = data.frame(org = tname, x = ti$x, fam = ti$cat3, score = tr[,tname])
for (qname in qnames) {
  tps = data.frame(org = qname, x = ti$x, fam = ti$cat3, score = tr[,qname])
  tp = rbind(tp, tps)
}

orgs = labels(r.dendro)
tp$org = factor(tp$org, levels = orgs)
p = ggplot(tp, aes(x = x, y = org, fill = score)) +
  geom_tile(stat = 'identity', position = "identity") + 
  scale_fill_gradientn(colours = matlab.like2(20)) +
#  scale_fill_manual(name = "AA distance:", breaks = labs, labels = labs, values = cols, guide = guide_legend(nrow = 2, byrow = T, label.position = "right", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
#  labs(fill = "Pariwise AA distance") +
#  coord_flip() +
  scale_x_discrete(name = '', breaks = tx$mid, labels = tx$cat) +
  scale_y_discrete(name = '') +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.1, 'lines')) +
  theme(legend.position = "right", legend.key.size = unit(1, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 8, colour = "brown", angle = 45, hjust = 0, vjust = 1)) +
  theme(axis.text.y = element_text(size = 10, colour = "blue", angle = 0))

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

fp = sprintf("%s/61.heatmap.%s.pdf", diro, fam)
pdf(file = fp, width = 15, height = 6, bg = 'transparent')
plot.new()

vl <- viewport(x = 0, y = 0, 
  width = unit(0.1, 'npc'), height = unit(1, 'npc'), 
  just = c('left', 'bottom'), name = 'left')
pushViewport(vl)
par(new = T, fig = gridFIG(), mar=c(4.3,0.4,0.2,0))
plot.phylo(as.phylo(r.hc), no.margin = F, show.tip.label = F)
upViewport()

vr <- viewport(x = unit(1, 'npc'), y = 0, 
  width = unit(0.9, 'npc'), height = unit(1, 'npc'), 
  just = c('right', 'bottom'), name = 'right')
pushViewport(vr)
grid.draw(gt)
upViewport()

dev.off()


##### plot NBS-LRR ortholog map
ft = file.path(diro, "01.ortho.tbl")
tt = read.table(ft, header = T, sep = "\t", as.is = T)

ti = tt[tt$cat == 'NBS-LRR',]
ti = merge(ti, tgene[,c(1,3)], by.x = tname, by.y = 'id')
ti = ti[order(ti$cat3, ti[,tname]),]
ti = cbind(ti, x = 1:nrow(ti))
tp = data.frame(org = tname, x = ti$x, fam = ti$cat3, refid = ti[,tname], score = 1)

for (qname in qnames) {
  fc = sprintf("%s/%s_%s/51_ortho/05.score.tbl", Sys.getenv("misc3"), qname, tname)
  tc = read.table(fc, sep = "\t", as.is = T, header = T)
  tc = cbind(tc[,c('tid','qid')], score = tc$mat/tc$len, stringsAsFactors = F)
  
  tis = ti[!is.na(ti[,qname]), c(tname, qname, 'x', 'cat3')]
  tps = merge(tis, tc, by.x = c(tname, qname), by.y = c('tid', 'qid'))
  tps = cbind(org = qname, tps[,c('x','cat3',tname,'score')])
  colnames(tps)[3:4] = c('fam', 'refid')
  tp = rbind(tp, tps)
}

tx = ddply(ti[,c('cat3','x')], .(cat3), summarise, bed = min(x), end = max(x), mid = as.integer(mean(x)))

tp$org = factor(tp$org, levels = rev(c(tname, qnames)))
p = ggplot(tp, aes(x = x, y = org, fill = score)) +
  geom_tile(stat = 'identity', position = "identity") + 
  scale_x_discrete(name = '', breaks = tx$mid, labels = tx$cat) +
  scale_y_discrete(name = '') +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 8, colour = "brown", angle = 45, hjust = 1, vjust = 1)) +
  theme(axis.text.y = element_text(size = 10, colour = "blue", angle = 0))

fp = sprintf("%s/nbs.pdf", diro)
ggsave(p, filename = fp, width = 15, height = 6)

##### plot misc GeneFam ortholog map
fams = c('2OG-FeII_Oxy', 'Peroxidase', 'Cytochrome', 'Hydrolase', 'Peptidase')
ft = file.path(diro, "01.ortho.tbl")
tt = read.table(ft, header = T, sep = "\t", as.is = T)

ti = tt[tt$cat %in% fams,]
ti = merge(ti, tgene[,c(1,3)], by.x = tname, by.y = 'id')
ti = ti[order(ti$cat3, ti[,tname]),]
ti = cbind(ti, x = 1:nrow(ti))
tp = data.frame(org = tname, x = ti$x, fam = ti$cat3, refid = ti[,tname], score = 1)

for (qname in qnames) {
  fc = sprintf("%s/%s_%s/51_ortho/05.score.tbl", Sys.getenv("misc3"), qname, tname)
  tc = read.table(fc, sep = "\t", as.is = T, header = T)
  tc = cbind(tc[,c('tid','qid')], score = tc$mat/tc$len, stringsAsFactors = F)
  
  tis = ti[!is.na(ti[,qname]), c(tname, qname, 'x', 'cat3')]
  tps = merge(tis, tc, by.x = c(tname, qname), by.y = c('tid', 'qid'))
  tps = cbind(org = qname, tps[,c('x','cat3',tname,'score')])
  colnames(tps)[3:4] = c('fam', 'refid')
  tp = rbind(tp, tps)
}

tx = ddply(ti[,c('cat','x')], .(cat), summarise, bed = min(x), end = max(x), mid = as.integer(mean(x)))

tp$org = factor(tp$org, levels = rev(c(tname, qnames)))
p = ggplot(tp, aes(x = x, y = org, fill = score)) +
  geom_tile(stat = 'identity', position = "identity") + 
  scale_x_discrete(name = '', breaks = tx$mid, labels = tx$cat) +
  scale_y_discrete(name = '') +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 8, colour = "brown", angle = 45, hjust = 1, vjust = 1)) +
  theme(axis.text.y = element_text(size = 10, colour = "blue", angle = 0))

fp = sprintf("%s/%s.pdf", diro, fams[1])
ggsave(p, filename = fp, width = 15, height = 6)

##### plot proteome diversity
qname = "HM058"
tname = "HM101"
dir = sprintf("%s/%s_%s/23_blat", Sys.getenv("misc3"), qname, tname)

fi = file.path(dir, "ortho.2.tbl")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

  tdir = sprintf("%s/%s", Sys.getenv("genome"), tname)
  fg = file.path(tdir, "51.gtb")
  tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1,16)]

tu = merge(ti, tg, by.x = 'tid', by.y = 'id', all.x = T)
colnames(tu)[7] = 'fam'

breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
tu = cbind(tu, dis = cut(tu$ident, breaks, include.lowest = T))
to = ddply(tu, .(fam, dis), summarise, cnt = length(fam))
to2 = ddply(tu, .(fam), summarise, cnt_fam = length(fam))
to = merge(to, to2, by = 'fam')
to = cbind(to, pct = to$cnt / to$cnt_fam)

tos = to[to$dis == '[0,0.01]',]
fams = tos$fam[order(tos$pct)]
to$fam = factor(to$fam, levels = fams)

x = as.numeric(to2$cnt_fam)
names(x) = to2$fam
cnts =  format(x[fams], big.mark = ",")

labs = unique(to$dis)
cols = rainbow(11)


p = ggplot(to, aes(x = fam, y = pct, fill = dis)) +
  geom_bar(stat = 'identity') + 
  scale_fill_manual(name = "AA distance:", breaks = labs, labels = labs, values = cols, guide = guide_legend(nrow = 2, byrow = T, label.position = "right", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
#  labs(fill = "Pariwise AA distance") +
  coord_flip() +
  scale_y_continuous(name = '', expand = c(0, 0)) +
  scale_x_discrete(name = '', expand = c(0, 0), labels = cnts) +
#  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0, 'lines')) +
  theme(legend.position = "bottom", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(2,2,0,5), "lines")) +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 8, colour = "grey", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0)) +
  annotation_custom(grob = textGrob(label = "Gene Family | #", just = c('right', 'bottom'), gp = gpar(fontsize = 8)), ymin = -0.05, ymax = -0.05, xmin = length(fams)+1, xmax = length(fams)+1)

for (i in 1:length(fams)) {
  p <- p + annotation_custom(grob = textGrob(label = fams[i], just = c('right'), gp = gpar(fontsize = 8)), ymin = -0.08, ymax = -0.08, xmin = i, xmax = i)
}
  
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

fp = sprintf("%s/compstat/%s.pdf", Sys.getenv("misc3"), tolower(qname))
pdf(file = fp, width = 8, height = 10, bg = 'transparent')
grid.newpage()
grid.draw(gt)
dev.off()

##### plot NBS-LRR / CRP diversity
qname = "HM004"
tname = "HM101"
dir = sprintf("%s/%s_%s/23_blat", Sys.getenv("misc3"), qname, tname)

fi = file.path(dir, "ortho.2.tbl")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

  tdir = sprintf("%s/%s", Sys.getenv("genome"), tname)
  fg = file.path(tdir, "51.gtb")
  tg = read.table(fg, sep = "\t", header = T, as.is = T)
  tgs = tg[tg$cat2 == 'CRP' | tg$cat2 == 'NBS-LRR', c('id','cat3')]

tu = merge(ti, tgs, by.x = 'tid', by.y = 'id', all.x = F)
colnames(tu)[7] = 'fam'

breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
tu = cbind(tu, dis = cut(tu$ident, breaks, include.lowest = T))
to = ddply(tu, .(fam, dis), summarise, cnt = length(fam))
to2 = ddply(tu, .(fam), summarise, cnt_fam = length(fam))
to = merge(to, to2, by = 'fam')
to = cbind(to, pct = to$cnt / to$cnt_fam)
to = to[to$cnt_fa >= 10,]

#tos = ddply(to, .(fam), summarise, pct_median = median(pct))
#fams = tos$fam[order(tos$pct_median)]
fams = sort(unique(to$fam))
to$fam = factor(to$fam, levels = fams)

x = as.numeric(to2$cnt_fam)
names(x) = to2$fam
cnts =  format(x[fams], big.mark = ",")

labs = unique(to$dis)
cols = rainbow(11)

p = ggplot(to, aes(x = fam, y = pct)) +
  geom_boxplot() + 
  coord_flip() +
  scale_y_continuous(name = '', expand = c(0, 0), limits = c(0,1)) +
  scale_x_discrete(name = '', expand = c(0, 0), labels = cnts) +
  theme(legend.position = "bottom", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(2,2,0,5), "lines")) +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 8, colour = "grey", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0)) +
  annotation_custom(grob = textGrob(label = "Gene Family | #", just = c('right', 'bottom'), gp = gpar(fontsize = 8)), ymin = -0.05, ymax = -0.05, xmin = length(fams)+1, xmax = length(fams)+1)

for (i in 1:length(fams)) {
  p <- p + annotation_custom(grob = textGrob(label = fams[i], just = c('right'), gp = gpar(fontsize = 8)), ymin = -0.08, ymax = -0.08, xmin = i, xmax = i)
}
  
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

fp = sprintf("%s/compstat/%s.crp+nbs.pdf", Sys.getenv("misc3"), tolower(qname))
pdf(file = fp, width = 8, height = 10, bg = 'transparent')
grid.newpage()
grid.draw(gt)
dev.off()
