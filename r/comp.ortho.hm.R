require(plyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(RColorBrewer)
require(reshape2)
source("Location.R")
source("comp.fun.R")

dirw = file.path(Sys.getenv("misc3"), "comp.ortho.hm")

## run comp.ortho.score.R to create 05_score/*.tbl
fg = file.path(tcfg$gdir, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:5,16)]

tg = tg[order(tg$chr, tg$beg, tg$end),]
tg = cbind(tg, idx = 1:nrow(tg))
colnames(tg)[6] = 'fam'
tg = rename_genefam(tg, fams)

## compute pairwise aln similarity (HM101-based)
to = tg[,1:2]
colnames(to)[1] = tname
for (qname in qnames_15) {
  fs = sprintf("%s/11_score/%s.tbl", dirw, qname)
  ts = read.table(fs, sep = "\t", header = T, as.is = T)
  ts2 = ts[(ts$qgap + ts$tgap) / ts$len <= 0.5,]
  ts3 = data.frame(tname = ts2$tid, qname = ts2$mat / (ts2$mat+ts2$mis), stringsAsFactors = F)
  colnames(ts3) = c(tname, qname)

  to = merge(to, ts3, by = tname, all = T)
#  system(sprintf("cp %s %s/01_syn_ortho/%s.tbl", fh, dirw, qname))
  cat(qname, nrow(to), "\n")
}
to = to[,-2]

norgs = apply(to, 1, function(x) sum(!is.na(x[-1])))
table(norgs)

fo = file.path(dirw, "12.score.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

##
fs = file.path(dirw, "12.score.tbl")
ts = read.table(fs, sep = "\t", header = T, as.is = T)
ids_nomissing = ts$HM101[norgs >= 10]

norgs = apply(ts, 1, function(x) sum(!is.na(x[-1])))
table(norgs)

## plot all pfam-misc
tgs = tg[tg$chr != 'chrU' & tg$id %in% ids_nomissing & !tg$fam %in% c("Unknown", "TE"),]
tgs = tgs[order(tgs$idx),]
tgs$idx = 1:nrow(tgs)

tx = ddply(tgs, .(chr), summarise, beg = min(idx), end = max(idx))
tx = tx[tx$end > tx$beg,]

tw = ts[ts$HM101 %in% tgs$id,]
colnames(tw)[1] = 'id'
tl = reshape(tw, direction = 'long', varying = list(2:16), idvar = 'id', timevar = 'org', v.names = 'score', times = colnames(tw)[2:16])

to = merge(tl, tgs[,c('id','idx')], by = 'id')
to$org = factor(to$org, levels = rev(qnames_15))

pb <- ggplot(to) +
  geom_tile(aes(x = idx, y = org, fill = 100*score), height = 0.8) +
  theme_bw() + 
  scale_x_continuous(name = '', limits = c(0, max(to$idx)+1), expand=c(0.01, 0), breaks = floor((tx$beg+tx$end)/2), labels = tx$chr) +
  scale_y_discrete(expand = c(0.02, 0), name = '') +
  scale_fill_gradient(name = 'sequence similarity with HM101 ortholog', space = "Lab", low = 'firebrick1', high = 'dodgerblue', na.value = 'grey50') +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0, 1), legend.title = element_text(size = 8), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,0.5,0,0), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid")) +
  annotate('segment', x = tx$beg, xend = tx$end, y = 0.3, yend = 0.3, size = 0.5, color = 'darkgreen') +
  annotate('segment', x = tx$beg, xend = tx$beg, y = 0.25, yend = 0.45, size = 0.5, color = 'darkgreen') +
  annotate('segment', x = tx$end, xend = tx$end, y = 0.25, yend = 0.45, size = 0.5, color = 'darkgreen')

fo = file.path(dirw, "13.score.pdf")
ggsave(pb, filename = fo, width = 8, height = 4)

## plot selected fams
famsp = c("Zinc-Finger", "CRP-NCR", "NBS-LRR")
pname = 'fams0'

plots = list()
for (i in 1:length(famsp)) {
  fam = famsp[i]

  tgs = tg[tg$chr != 'chrU' & tg$id %in% ids_nomissing & tg$fam == fam,]
  tgs = tgs[order(tgs$idx),]
  tgs$idx = 1:nrow(tgs)

  tx = ddply(tgs, .(chr), summarise, beg = min(idx), end = max(idx))
  tx = tx[tx$end > tx$beg,]

  tw = ts[ts$HM101 %in% tgs$id,]
  colnames(tw)[1] = 'id'
  tl = reshape(tw, direction = 'long', varying = list(2:16), idvar = 'id', timevar = 'org', v.names = 'score', times = colnames(tw)[2:16])

  to = merge(tl, tgs[,c('id','idx')], by = 'id')
  to$org = factor(to$org, levels = rev(qnames_15))

p1 <- ggplot(to) +
  geom_tile(aes(x = idx, y = org, fill = 100*score), height = 0.8) +
  theme_bw() + 
  scale_x_continuous(name = '', limits = c(0, max(to$idx)+1), expand=c(0.01, 0), breaks = floor((tx$beg+tx$end)/2), labels = tx$chr) +
  scale_y_discrete(expand = c(0.02, 0), name = '') +
  scale_fill_gradient(name = 'sequence similarity with HM101 ortholog', space = "Lab", low = 'firebrick1', high = 'dodgerblue', na.value = "white") +
  theme(plot.margin = unit(c(1,0.5,0,0), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid")) +
  annotate('segment', x = tx$beg, xend = tx$end, y = 0.3, yend = 0.3, size = 0.5, color = 'darkgreen') +
  annotate('segment', x = tx$beg, xend = tx$beg, y = 0.25, yend = 0.45, size = 0.5, color = 'darkgreen') +
  annotate('segment', x = tx$end, xend = tx$end, y = 0.25, yend = 0.45, size = 0.5, color = 'darkgreen')

  if(i == 1) {
    p1 <- p1 + theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0, 1), legend.title = element_text(size = 8), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill=NA, size=0), legend.margin = unit(0, "line"))
  } else {
    p1 <- p1 + theme(legend.position = 'none')
  }
  plots[[fam]] = p1
}

numrow = length(famsp)

fo = sprintf("%s/13.score.%s.pdf", dirw, pname)
pdf(file = fo, width = 6, height = 2*numrow+0.5, bg = 'transparent')
grid.newpage()
pushViewport(viewport(layout = grid.layout(numrow, 1, height = c(2.5,rep(2, numrow-1)))))

dco = data.frame(x = 1:numrow, y = rep(1,numrow), lab = paste(LETTERS[1:numrow], famsp, sep = ": "))
for (i in 1:nrow(dco)) {
  x = dco$x[i]; y = dco$y[i]; lab = dco$lab[i]
  print(plots[[famsp[i]]], vp = viewport(layout.pos.row = x, layout.pos.col = y))
  grid.text(lab, x = 0, y = unit(1, 'npc'), just = c('left', 'top'), 
    gp =  gpar(col = "black", fontface = 2, fontsize = 15),
    vp = viewport(layout.pos.row = x, layout.pos.col = y))
}
dev.off()


