require(plyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(RColorBrewer)
source('comp.fun.R')

dirw = file.path(Sys.getenv("misc3"), "comp.ortho")
diro = file.path(Sys.getenv("misc3"), "comp.stat")

fi = file.path(dirw, "21.gid.tbl")
tid = read.table(fi, header = T, sep = "\t", as.is = T)
fi = file.path(dirw, "21.loc.tbl")
tlo = read.table(fi, header = T, sep = "\t", as.is = T)
fi = file.path(dirw, "21.sta.tbl")
tst = read.table(fi, header = T, sep = "\t", as.is = T)
fi = file.path(dirw, "22.cat.tbl")
tca = read.table(fi, header = T, sep = "\t", as.is = T)
fi = file.path(dirw, "28.dist.tbl")
tds = read.table(fi, header = T, sep = "\t", as.is = T)

norgs = apply(ti[,orgs], 1, function(x) sum(!x %in% c('','-')))
table(norgs)

##### gene fam mpd distribution 
orgs_in = c(tname, get_orgs("ingroup"))
  comps = c()
  ids = orgs_in
  for (i in 1:(length(ids)-1)) {
    id1 = ids[i]
    for (j in (i+1):length(ids)) {
      id2 = ids[j]
      comps = c(comps, sprintf("%s.%s", id1, id2))
    }
  }

norgs = apply(tst[,orgs_in], 1, function(x) sum(x %in% c('syn','rbh')))
idxs = which(norgs > 1)
mpd = apply(tds[idxs,comps], 1, mean, na.rm = T)
idxs = idxs[!is.na(mpd)]
mpd = mpd[!is.na(mpd)]
tu = data.frame(fam = tca$cat2[idxs], mpd = mpd, stringsAsFactors = F)

breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
tu = cbind(tu, dis = cut(tu$mpd, breaks, include.lowest = T))
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
  theme(legend.position = "bottom", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(2,1,0,5), "lines")) +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 8, colour = "grey", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0)) +
  annotation_custom(grob = textGrob(label = "Gene Family | #", just = c('right', 'bottom'), gp = gpar(fontsize = 8)), ymin = -0.05, ymax = -0.05, xmin = length(fams)+1, xmax = length(fams)+1)

for (i in 1:length(fams)) {
  p <- p + annotation_custom(grob = textGrob(label = fams[i], just = c('right'), gp = gpar(fontsize = 8)), ymin = -0.08, ymax = -0.08, xmin = i, xmax = i)
}
  
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

fp = file.path(diro, "76.proteome.mpd.pdf")
pdf(file = fp, width = 8, height = 10, bg = 'transparent')
grid.newpage()
grid.draw(gt)
dev.off()


##### gene fam AFS
norgs = apply(tst[,orgs], 1, function(x) sum(x %in% c('syn','rbh')))
fam_cnt = table(tca$cat2)
to0 = cbind(tca, n_org = norgs, stringsAsFactors = F)
to1 = ddply(to0, .(n_org, cat2), summarise, cnt = length(cat2))
to2 = cbind(to1, fam_cnt = fam_cnt[to1$cat2])
to = cbind(to2, fam_pct = to2$cnt / to2$fam_cnt)

tos = to[to$n_org == 16,]
fams = tos$cat2[order(tos$fam_pct)]
all_fams = unique(to$cat2)
fams = c(fams, all_fams[!all_fams %in% fams])
to$cat2 = factor(to$cat2, levels = fams)
cnts = format(fam_cnt[fams], big.mark = ',')

p = ggplot(to, aes(x = cat2, y = fam_pct, fill = n_org)) +
  geom_bar(stat = 'identity', position = "stack") + 
#  scale_fill_manual(name = "AA distance:", breaks = labs, labels = labs, values = cols, guide = guide_legend(nrow = 2, byrow = T, label.position = "right", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  labs(fill = "Sharing Accessions") +
  coord_flip() +
  scale_x_discrete(name = '', labels = cnts) +
  scale_y_continuous(name = '', expand = c(0, 0)) +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.1, 'lines')) +
  theme(legend.position = "top", legend.key.size = unit(0.6, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,1,0,4.5), "lines")) +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1)) +
  annotation_custom(grob = textGrob(label = "Gene Family | #", just = c('right', 'bottom'), gp = gpar(fontsize = 8)), ymin = -0.04, ymax = -0.04, xmin = length(fams)+1, xmax = length(fams)+1)

for (i in 1:length(fams)) {
  p <- p + annotation_custom(grob = textGrob(label = fams[i], just = c('right'), gp = gpar(fontsize = 8)), ymin = -0.07, ymax = -0.07, xmin = i, xmax = i)
}

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

fp = sprintf("%s/76.proteome.afs.pdf", diro)
pdf(file = fp, width = 7, height = 10, bg = 'transparent')
grid.newpage()
grid.draw(gt)
dev.off()

##### plot selected fams / NBS-LRR / CRP AFS
fams = c("Peroxidase", "Hydrolase", "Auxin_inducible", "Deaminase",
  "tnl0850", "tnl0480", "tnl0400", "tnl0100",
  "cnl1500", "cnl1400", "cnl0400", "cnl0200", "cnl0950", 
  "CRP0010", "CRP0110", "CRP0355", "CRP0675", 
  "CRP1190", "CRP1230", "CRP1430", "CRP1520", "CRP1530")
fams = sortC(unique(fams))
norgs = apply(tst[,orgs], 1, function(x) sum(x %in% c('syn','rbh')))
idxs = which(tca$cat3 %in% fams | tca$cat2 %in% fams)

tx = data.frame(idx = idxs, fam = tca$cat2[idxs], cat3 = tca$cat3[idxs], norg = norgs[idxs], stringsAsFactors = F)
idxs2 = which(tx$fam %in% c("CRP", "NBS-LRR"))
tx$fam[idxs2] = tx$cat3[idxs2]

tx1 = ddply(tx, .(fam, norg), summarise, cnt = length(fam))
tx2 = ddply(tx, .(fam), summarise, cnt_fam = length(fam))
to = merge(tx1, tx2, by = 'fam')
to = cbind(to, pct = to$cnt / to$cnt_fam)
to = to[to$cnt_fam >= 10,]

th = ddply(to, .(fam), summarise, het = 1 - sum(pct^2))

fam_cnt = tx2$cnt_fam; names(fam_cnt) = tx2$fam
to$fam = factor(to$fam, levels = fams)
cnts = format(fam_cnt[fams], big.mark = ',')

p = ggplot(to, aes(x = fam, y = pct, fill = norg)) +
  geom_bar(stat = 'identity', position = "stack") + 
#  scale_fill_manual(name = "AA distance:", breaks = labs, labels = labs, values = cols, guide = guide_legend(nrow = 2, byrow = T, label.position = "right", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  labs(fill = "Sharing Accessions") +
  coord_flip() +
  scale_x_discrete(name = '', labels = cnts) +
  scale_y_continuous(name = '', expand = c(0, 0)) +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.1, 'lines')) +
  theme(legend.position = "top", legend.key.size = unit(0.6, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,1,0,4.5), "lines")) +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1)) +
  annotation_custom(grob = textGrob(label = "Gene Family | #", just = c('right', 'bottom'), gp = gpar(fontsize = 8)), ymin = -0.03, ymax = -0.03, xmin = length(fams)+1, xmax = length(fams)+1)

for (i in 1:length(fams)) {
  p <- p + annotation_custom(grob = textGrob(label = fams[i], just = c('right'), gp = gpar(fontsize = 8)), ymin = -0.07, ymax = -0.07, xmin = i, xmax = i)
}

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

fp = sprintf("%s/77.afs.select.pdf", diro)
pdf(file = fp, width = 6, height = 4, bg = 'transparent')
grid.newpage()
grid.draw(gt)
dev.off()
