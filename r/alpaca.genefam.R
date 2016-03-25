require(plyr)
require(dplyr)
require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(xlsx)
require(seqinr)
require(RColorBrewer)
require(ape)
source("Location.R")
source("comp.fun.R")

diro = file.path(Sys.getenv('misc3'), 'alpaca')
setwd(diro)

qnames = qnames_alpaca_comp

tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

## create a config file
fps = sprintf("%s/%s/51.fas", Sys.getenv('genome'), qnames)
to = data.frame(org = qnames, fp = fps, stringsAsFactors = F)
fo = file.path(diro, "31.conf.csv")
#write.table(to, fo, sep = ",", row.names = F, col.names = F, quote = F)
#merge.fas.py 31.conf.csv 32.pro.fas

## read all genes
tg = data.frame()
for (qname in qnames) {
  dirg = sprintf("%s/%s", Sys.getenv("genome"), qname)
  fg = file.path(dirg, "51.gtb")
  tg1 = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1,3:6,16:17)]
  tg = rbind(tg, cbind(org = qname, tg1))
}
colnames(tg)[c(2,7:8)] = c('gid','fam','sfam')

## write gene ids 
sfams = c('CRP0355', 'CRP3710', 'CRP4180')
opts = c('alps', 'alpc')

for (sfam in sfams) {
  tgs = tg[tg$sfam == sfam,]
  tgs = cbind(tgs, id = sprintf("%s|%s", tgs$org, tgs$gid))
  ids1 = tgs$id[tgs$org %in% c("HM101", "HM056", "HM034", "HM340")]
  ids2 = tgs$id[tgs$org %in% c("HM101", "HM056.AC", "HM034.AC", "HM340.AC")]
  fo1 = sprintf("%s/34.sfams/%s.alps.bed", diro, sfam)
  fo2 = sprintf("%s/34.sfams/%s.alpc.bed", diro, sfam)
  write.table(ids1, fo1, sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(ids2, fo2, sep = "\t", row.names = F, col.names = F, quote = F)
}

### run alpaca.genefam.py
cols_raw = brewer.pal(10, "Paired")
cols_map = cols_raw[c(1:6,10)]
names(cols_map) = qnames

f_cl = file.path(diro, "34.color.legend.pdf")
pdf(f_cl, width = 2, height = 2, bg = 'transparent')
par(mar=c(0,0,0,0))
plot(c(0, 200), c(0, 200), type= "n", axes = F)
legend(0, 200, legend = qnames, fill = cols_map, bty = 'n')
dev.off()

sfam = sfams[1]
for (opt in opts) {
  ft1 = sprintf("%s/34.sfams/%s.%s.phy.nwk", diro, sfam, opt)
  tree = read.tree(ft1)
  fo1 = sprintf("%s/34.sfams/%s.%s.png", diro, sfam, opt)
  
  res = strsplit(tree$tip.label, '[|]')
  oids = sapply(res, "[", 1)
  chrs = sapply(strsplit(sapply(res, "[", 2), "_"), "[", 2)
  tree$tip.label = chrs
  font = rep(1, length(oids))
  
  tip.cols = cols_map[oids]

  scores = as.numeric(tree$node.label)
  node.bg = rep('white', tree$Nnode)
  node.bg[scores >= 0.95] = 'black'
  node.bg[scores >= 0.8 & scores < 0.95] = 'gray'
  
  png(filename = fo1, width = 600, height = 800, units = 'px')
  plot(tree, type = 'phylogram', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols, label.offset = 0.08, font = font,
    no.margin = T, cex = 0.8)
  tiplabels(pch = 22, frame = 'none', adj = 0.55, bg = tip.cols)
 # nodelabels(pch = 22, bg = node.bg)
  add.scale.bar(x = 0.02, y = tree$Nnode*0.9 , lcol = 'black')
  dev.off()
}

sfam = sfams[2]
for (opt in opts) {
  ft1 = sprintf("%s/34.sfams/%s.%s.phy.nwk", diro, sfam, opt)
  tree = read.tree(ft1)
  fo1 = sprintf("%s/34.sfams/%s.%s.png", diro, sfam, opt)
  
  res = strsplit(tree$tip.label, '[|]')
  oids = sapply(res, "[", 1)
  chrs = sapply(strsplit(sapply(res, "[", 2), "_"), "[", 2)
  tree$tip.label = chrs
  font = rep(1, length(oids))
  
  tip.cols = cols_map[oids]

  scores = as.numeric(tree$node.label)
  node.bg = rep('white', tree$Nnode)
  node.bg[scores >= 0.95] = 'black'
  node.bg[scores >= 0.8 & scores < 0.95] = 'gray'
  
  png(filename = fo1, width = 600, height = 800, units = 'px')
  plot(tree, type = 'phylogram', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols, label.offset = 0.08, font = font,
    no.margin = T, cex = 0.8)
  tiplabels(pch = 22, frame = 'none', adj = 0.55, bg = tip.cols)
 # nodelabels(pch = 22, bg = node.bg)
  add.scale.bar(x = 0.02, y = tree$Nnode*0.9 , lcol = 'black')
  dev.off()
}

sfam = sfams[3]
for (opt in opts) {
  ft1 = sprintf("%s/34.sfams/%s.%s.phy.nwk", diro, sfam, opt)
  tree = read.tree(ft1)
  fo1 = sprintf("%s/34.sfams/%s.%s.png", diro, sfam, opt)
  
  res = strsplit(tree$tip.label, '[|]')
  oids = sapply(res, "[", 1)
  chrs = sapply(strsplit(sapply(res, "[", 2), "_"), "[", 2)
  tree$tip.label = chrs
  font = rep(1, length(oids))
  
  tip.cols = cols_map[oids]

  scores = as.numeric(tree$node.label)
  node.bg = rep('white', tree$Nnode)
  node.bg[scores >= 0.95] = 'black'
  node.bg[scores >= 0.8 & scores < 0.95] = 'gray'
  
  png(filename = fo1, width = 600, height = 800, units = 'px')
  plot(tree, type = 'phylogram', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols, label.offset = 0.06, font = font,
    no.margin = T, cex = 1)
  tiplabels(pch = 22, frame = 'none', adj = 0.52, bg = tip.cols)
 # nodelabels(pch = 22, bg = node.bg)
  add.scale.bar(x = 0.02, y = tree$Nnode*0.9 , lcol = 'black')
  dev.off()
}