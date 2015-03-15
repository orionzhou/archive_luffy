require(rtracklayer)
require(plyr)
require(grid)
require(ggplot2)
source('Location.R')
source('comp.fun.R')

dirw = sprintf("%s/comp.stat", Sys.getenv("misc3"))

tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

grnp = GenomicRanges::setdiff(grt, grp)
tsize = sum(width(grt))
tsize2 = sum(width(grnp))

tg = read.table(tcfg$gene, sep = "\t", header = F, as.is = T)
colnames(tg) = c('chr', 'beg', 'end', 'srd', 'id', 'type', 'cat')
tg = tg[tg$type == 'cds',]
#tg_size = ddply(tg, .(id, cat), summarise, size = sum(end - beg + 1))
grgt = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

##### SV impact on gene families
qname = "HM004"
for (qname in qnames_all) {

cfg = cfgs[[qname]]

tlen = read.table(cfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))

tgap = read.table(cfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

grnp = GenomicRanges::setdiff(grt, grp)
qsize = sum(width(grt))
qsize2 = sum(width(grnp))

qg = read.table(cfg$gene, sep = "\t", header = F, as.is = T)
colnames(qg) = c('chr', 'beg', 'end', 'srd', 'id', 'type', 'cat')
qg = qg[qg$type == 'cds',]
#qg_size = ddply(qg, .(id, cat), summarise, size = sum(end - beg + 1))
grgq = with(qg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

cdir = cfg$cdir
fy = file.path(cdir, "31.9/gal")
ty = read.table(fy, header = T, sep = "\t", as.is = T)[,c('tId', 'tBeg', 'tEnd', 'qId', 'qBeg', 'qEnd')]
gryt = with(ty, GRanges(seqnames = tId, ranges = IRanges(tBeg, end = tEnd)))
gryq = with(ty, GRanges(seqnames = qId, ranges = IRanges(qBeg, end = qEnd)))
ytlen = sum(width(reduce(gryt)))
yqlen = sum(width(reduce(gryq)))

fv = file.path(cdir, "../31_sv/05.stb")
tv = read.table(fv, header = T, sep = "\t", as.is = T)
tv = tv[,c('id','tchr','tbeg','tend','tlen','srd','qchr','qbeg','qend','qlen')]
tv = tv[tv$tlen+tv$qlen-2 >= 50,]
tvt = tv[tv$tlen > 0, c('tchr','tbeg','tend')]
grst = with(tvt, GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend)))
tvq = tv[tv$qlen > 0, c('qchr','qbeg','qend')]
grsq = with(tvq, GRanges(seqnames = qchr, ranges = IRanges(qbeg, end = qend)))

stlen = sum(width(reduce(grst)))
sqlen = sum(width(reduce(grsq)))

tocnt = intersect_count(grgt, grst)
tg2 = cbind(tg, cnt = tocnt)
tgb = group_by(tg2, id)
x = summarise(tgb, cnt = sum(cnt))
tpct = sum(x$cnt > 0) / nrow(x)

qocnt = intersect_count(grgq, grsq)
qg2 = cbind(qg, cnt = qocnt)
qgb = group_by(qg2, id)
x = summarise(qgb, cnt = sum(cnt))
qpct = sum(x$cnt > 0) / nrow(x)

cat(sprintf("%s: genome [tgt %.03f qry %.03f] gene [tgt %.03f qry %.03f]\n", qname, stlen/ytlen, sqlen/yqlen, tpct, qpct))

}

fo = sprintf("%s/comp_genesv_%s.pdf", dirw)
p = ggplot(do) +
  geom_bar(mapping = aes(x = cat, y = mean), 
    stat = 'identity', geom_params=list(width = 0.8, alpha = 0.8)) +
  geom_errorbar(mapping = aes(x = cat, ymin = mean-std, ymax = mean+std), 
    stat = 'identity', geom_params=list(width = 0.4, alpha = 0.8)) +
#  scale_fill_brewer(palette='Set1', name = '', guide = guide_legend(nrow = 1, byrow = T, label.position = "right", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  coord_flip() +
  scale_x_discrete(name = 'fam') +
  scale_y_continuous(name = 'pct') +
  facet_wrap( ~ type, ncol = 4) +
  theme_bw() +
  theme(legend.position = "top", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(0,1,1,0), "lines")) +
  theme(axis.text.x = element_text(size = 7, colour = "gray", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0))
ggsave(p, filename = fo, width = 8, height = 10)


