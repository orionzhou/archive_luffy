require(plyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(RColorBrewer)

tname = "HM101"
qnames = c(
  "HM058", "HM125", "HM056.AC", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340.AC"
)
diro = file.path(Sys.getenv("misc3"), "comp.ortho")

fi = file.path(diro, "33.ortho.cat.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti[is.na(ti)] = ''
ti = ti[ti$cat2 != 'TE',]

orgs = c(tname, qnames)
n_org = apply(ti, 1, function(z) sum(z[1:length(orgs)] != ''))
ti = cbind(ti, n_org = n_org)

##### plot pan-proteome AFS
tab1 = table(ti$n_org)
tab1 = tab1[names(tab1) != 1]
dt1 = data.frame(n_org = as.numeric(names(tab1)), cnt = as.numeric(tab1), org = 'mixed', stringsAsFactors = F)

tis = ti[n_org == 1,]
oorg = apply(tis, 1, function(z) orgs[which(z[1:length(orgs)] != '')])
tab2 = table(oorg)
dt2 = data.frame(n_org = 1, cnt = as.numeric(tab2), org = names(tab2), stringsAsFactors = F)

dt = rbind(dt1, dt2)

cols = c(brewer.pal(12, 'Set3'), brewer.pal(4, 'Set1'), 'gray30')
labs = orgs

dt$org = factor(dt$org, levels = c(orgs, 'mixed'))
dt$n_org = factor(dt$n_org, levels = sort(as.numeric(unique(dt$n_org))))
p = ggplot(dt, aes(x = n_org, y = cnt, fill = org, order = plyr:::desc(org))) +
  geom_bar(stat = 'identity', position = "stack", geom_params=list(width = 0.5)) +
  scale_fill_manual(name = "Accession-Specific", breaks = labs, labels = labs, values = cols, guide = guide_legend(ncol = 1, byrow = F, label.position = "right", direction = "vertical", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  scale_x_discrete(name = '# Sharing Accession') +
  scale_y_continuous(name = '', expand = c(0.02, 0)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(legend.position = "right", legend.key.size = unit(0.6, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 7, colour = "blue")) +
  theme(axis.text.y = element_text(size = 7, colour = "brown", angle = 0, hjust = 1))

fp = sprintf("%s/71.ortho.afs.pdf", diro)
ggsave(p, filename = fp, width = 7, height = 5)


##### pan-proteome size
reps = 1:10
n_orgs = 1:16
tp = data.frame(rep = rep(reps, each = length(n_orgs)), 
  n_org = rep(n_orgs, length(rep)), core = NA, pan = NA)

for (rep in reps) {
  set.seed(rep * 100)
  orgs = c(tname, sample(qnames))
  for (i in 1:length(orgs)) {
    orgss = orgs[1:i]
    tis = data.frame(ti[, orgss])
    pan = sum(apply(tis, 1, function(z) sum(z != '') >= 1))
    core = sum(apply(tis, 1, function(z) sum(z == '') == 0))
    tp$pan[tp$rep == rep & tp$n_org == i] = pan
    tp$core[tp$rep == rep & tp$n_org == i] = core
    cat(rep, orgs[i], core, pan, "\n")
  }
}

tp$rep = factor(tp$rep, levels = 1:max(tp$rep))
p = ggplot(tp) +
  geom_point(aes(x = n_org, y = pan), shape = 1, size = 1.1) +
  geom_point(aes(x = n_org, y = core), shape = 4, size = 1.1) +
  stat_smooth(aes(x = n_org, y = pan, col = 'pan'), fill = 'azure4', size = 0.2) +
  stat_smooth(aes(x = n_org, y = core, col = 'core'), fill = 'azure4', size = 0.2) +
#  scale_shape(name = "", solid = FALSE, guide = F) +
  scale_color_manual(name = "", labels = c('Core-proteome', 'Pan-proteome'), values = c("dodgerblue", "firebrick1"), guide = guide_legend(label.position = "left", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  scale_x_continuous(name = '# Genomes Sequenced') +
  scale_y_continuous(name = '# Genes', expand = c(0.02, 0), limits = c(0, 180000)) + 
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.2, 'lines')) +
  theme(legend.position = "top", legend.key.size = unit(1, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,1,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, color = "blue")) +
  theme(axis.text.y = element_text(size = 8, color = "grey", angle = 90, hjust  = 0.5))
  
fp = sprintf("%s/73.pan.size.pdf", diro)
ggsave(p, filename = fp, width = 6, height = 5)

##### plot genefam AFS
fam_cnt = table(ti$cat2)
to1 = ddply(ti, .(n_org, cat2), summarise, cnt = length(cat2))
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
  theme(legend.position = "top", legend.key.size = unit(0.6, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(1,1,0,4.5), "lines")) +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))

for (i in 1:length(fams)) {
  p <- p + annotation_custom(grob = textGrob(label = fams[i], just = c('right'), gp = gpar(fontsize = 8)), ymin = -0.07, ymax = -0.07, xmin = i, xmax = i)
}

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

fp = sprintf("%s/81.genefam.afs.pdf", diro)
pdf(file = fp, width = 7, height = 10, bg = 'transparent')
grid.newpage()
grid.draw(gt)
dev.off()
