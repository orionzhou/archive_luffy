require(plyr)
require(ggplot2)
require(seqinr)
require(GenomicRanges)

diro = file.path(Sys.getenv("misc3"), "comp.ortho")

tname = "HM101"
qnames = c(
  "HM058", "HM125", "HM056", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340"
)

f_tgene = file.path(Sys.getenv("genome"), tname, "51.gtb")
tgene = read.table(f_tgene, header = T, sep = "\t", as.is = T)[,c(1,16:17)]

##### create ortholog groups for 16 accessions
to = tgene[,1:2]
colnames(to) = c(tname, 'cat')
i = 1
for (qname in qnames) {
  fi = sprintf("%s/%s_%s/51_ortho/31.ortho.tbl", Sys.getenv("misc3"), qname, tname)
  ti = read.table(fi, sep = "\t", as.is = T, header = T)
  to = merge(to, ti, by.x = tname, by.y = 'tid', all.x = T)
  colnames(to)[2+i] = qname
  i = i + 1
}

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
dir = sprintf("%s/%s_%s/23_blat", Sys.getenv("misc3"), toupper(qname), toupper(tname))

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
