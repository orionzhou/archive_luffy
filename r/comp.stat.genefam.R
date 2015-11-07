require(plyr)
require(ggplot2)
require(GenomicRanges)
require(grid)
require(RColorBrewer)
require(gtable)
source('comp.fun.R')

dirw = file.path(Sys.getenv("misc3"), "comp.ortho")
diro = file.path(Sys.getenv("misc3"), "comp.stat")

ffam = "/home/youngn/zhoup/Data/db/pfam/genefam.tbl"
tfam = read.table(ffam, header = T, sep = "\t", as.is = T)[,1:2]

fams = c("CC-NBS-LRR", "TIR-NBS-LRR", "F-box", "LRR-RLK", "NCR", "TE", "Unknown")
fams = c(fams, "CRP0000-1030", "CRP1600-6250", "RLK")

tg = read.table(tcfg$gene, sep = "\t", header = F, as.is = T)
colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
tg = tg[tg$type == 'cds',]

##### genefam three-panel plot: thetaPi + largeEff + cdslen + mpd + cnv
# pi
ff = file.path(diro, "42.pi.tbl")
tf = read.table(ff, header = T, sep = "\t", as.is = T)

to = tf
to$fam[! to$fam %in% fams] = 'Pfam-Miscellaneous'
to$fam = factor(to$fam, levels = c(fams, 'Pfam-Miscellaneous'))

tp = ddply(to, .(fam), summarise, cnt = length(fam), q25 = quantile(pi, 0.25, na.rm = T), q50 = median(pi, na.rm = T), q75 = quantile(pi, 0.75, na.rm = T))
famsf = as.character(tp$fam[order(tp$q50, decreasing = T)])
tp$fam = factor(tp$fam, levels = famsf)
labs = sprintf("%s (%d)", famsf, tp$cnt[match(famsf, tp$fam)])

p1 = ggplot(tp) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', position = 'dodge', geom_params = list(width = 0.7)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01), labels = labs) +
  scale_y_continuous(name = 'ThetaPi', expand = c(0, 0), limits = c(0, 0.02)) +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

# largeeff
fi = file.path(diro, "43.largeeff.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)

to = ti
to$fam[! to$fam %in% fams] = 'Pfam-Miscellaneous'
to$fam = factor(to$fam, levels = c(fams, 'Pfam-Miscellaneous'))

tp = ddply(to, .(fam, eff), summarise, prop = sum(cnt > 0) / length(cnt))
tp$fam = factor(tp$fam, levels = famsf)

p2 = ggplot(tp) +
  geom_bar(aes(x = fam, y = prop, fill = eff),
    stat = 'identity', position = 'stack', geom_params = list(width = 0.7)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01)) +
  scale_y_continuous(name = 'Proportion', expand = c(0, 0), limits = c(0, 1)) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  theme(legend.position = c(0.7, 0.8), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.7, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_blank(), legend.text = element_text(size = 8, angle = 0)) +
  theme(plot.margin = unit(c(1,1,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
#  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))
  theme(axis.text.y = element_blank())

# protein length
gb = dplyr::group_by(tg, id)
to = dplyr::summarise(gb, fam = fam[1], len = sum(end-beg+1))

to$fam[! to$fam %in% fams] = 'Pfam-Miscellaneous'
to$fam = factor(to$fam, levels = c(fams, 'Pfam-Miscellaneous'))

tp = ddply(to, .(fam), summarise, q25 = quantile(len, 0.25), q50 = quantile(len, 0.5), q75 = quantile(len, 0.75))
tp$fam = factor(tp$fam, levels = famsf)
labs = sprintf("%s (%d)", famsf, tp$cnt[match(famsf, tp$fam)])

p3 = ggplot(tp) +
  geom_crossbar(aes(x = fam, y = q50/3, ymin = q25/3, ymax = q75/3),
    stat = 'identity', position = 'dodge', geom_params = list(width = 0.6)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01)) +
  scale_y_continuous(name = 'Protein Length', expand = c(0, 0), limits = c(0, 1300)) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  theme(plot.margin = unit(c(1,1,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_blank())

# mpd
fi = file.path(diro, "46.mpd.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)

to = ti
to$fam[! to$fam %in% fams] = 'Pfam-Miscellaneous'
to$fam = factor(to$fam, levels = c(fams, 'Pfam-Miscellaneous'))

tp = ddply(to, .(fam), summarise, cnt = length(fam), q25 = quantile(mpd, 0.25), q50 = quantile(mpd, 0.5), q75 = quantile(mpd, 0.75))
tp$fam = factor(tp$fam, levels = famsf)
labs = sprintf("%d", tp$cnt[match(famsf, tp$fam)])

p4 = ggplot(tp) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', position = 'dodge', geom_params = list(width = 0.6)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01), labels = labs) +
  scale_y_continuous(name = 'Mean Pariwise Protein Distance', expand = c(0,0), limits = c(0, 0.2)) +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

# cnv
fi = file.path(diro, "48.cnv.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)

to = ti
to$fam[! to$fam %in% fams] = 'Pfam-Miscellaneous'
to$fam = factor(to$fam, levels = c(fams, 'Pfam-Miscellaneous'))

tp = ddply(to, .(fam), summarise, cnt = length(fam), q25 = quantile(cv, 0.25), q50 = quantile(cv, 0.5), q75 = quantile(cv, 0.75))
tp$fam = factor(tp$fam, levels = famsf)
labs = sprintf("%d", tp$cnt[match(famsf, tp$fam)])
tp = rbind(tp, data.frame(fam = 'TE', cnt = 0, q25 = NA, q50 = NA, q75 = NA, stringsAsFactors = F))

p5 = ggplot(tp) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', geom_params = list(width = 0.5)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01), labels = labs) +
  scale_y_continuous(name = 'C.V. of cluster size') +
  theme_bw() +
#  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  theme(plot.margin = unit(c(1,1,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))
  

fo = file.path(diro, "49.genefam.pdf")
numcol = 5
wds = c(4, 3, 2.5, 3,3)
pdf(file = fo, width = sum(wds), height = 3.2, bg = 'transparent')
#tiff(file = fo, width = sum(wds), height = 4, units = 'in', bg = 'white')
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, numcol, width = wds)))

print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(p4, vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
print(p5, vp = viewport(layout.pos.row = 1, layout.pos.col = 5))

dco = data.frame(x = rep(1,numcol), y = 1:numcol, lab = LETTERS[1:numcol])
for (i in 1:nrow(dco)) {
  x = dco$x[i]; y = dco$y[i]; lab = dco$lab[i]
  grid.text(lab, x = 0, y = unit(1, 'npc'), just = c('left', 'top'), gp = gpar(col = "black", fontface = 2, fontsize = 20),
    vp = viewport(layout.pos.row = x, layout.pos.col = y))
}
dev.off()


#####
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

norgs = apply(tst[,qnames_15], 1, function(x) sum(!x %in% c('','-')))
table(norgs)

##### gene fam mpd distribution 
orgs_in = c(tname, qnames_12)
  comps = c()
  ids = orgs_in
  for (i in 1:(length(ids)-1)) {
    id1 = ids[i]
    for (j in (i+1):length(ids)) {
      id2 = ids[j]
      comps = c(comps, sprintf("%s.%s", id1, id2))
    }
  }

#### cross-bar plot showing quantile
norgs = apply(tst[,qnames_15], 1, function(x) sum(!x %in% c('','-')))
idxs = which(norgs > 10)
mpd = apply(tds[idxs,comps], 1, mean, na.rm = T)
idxs = idxs[!is.na(mpd)]
mpd = mpd[!is.na(mpd)]
tu = data.frame(fam = tca$cat2[idxs], sfam = tca$cat3[idxs], mpd = mpd, stringsAsFactors = F)

fo = file.path(diro, "46.mpd.tbl")
write.table(tu, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

tis = ddply(tu, .(fam), summarise, cnt = length(fam), q25 = quantile(mpd, 0.25), q50 = quantile(mpd, 0.5), q75 = quantile(mpd, 0.75))
tis = tis[order(tis$q50, decreasing = T),]
fams = tis$fam

tis$fam = factor(tis$fam, levels = fams)
p1 = ggplot(tis) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', position = 'dodge', geom_params = list(width = 0.6)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01), labels = sprintf("%s | %5d", tis$fam, tis$cnt)) +
  scale_y_continuous(name = 'Mean Pariwise Protein Distance', expand = c(0,0), limits = c(0, 0.2)) +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = file.path(diro, "46.mpd.pdf")
ggsave(p1, filename = fp, width = 5, height = 8)


## heatmap adopted from Arabidopsis paper (obsolete)
norgs = apply(tst[,orgs], 1, function(x) sum(!x %in% c('','-')))
idxs = which(norgs > 6)
mpd = apply(tds[idxs,comps], 1, mean, na.rm = T)
idxs = idxs[!is.na(mpd)]
mpd = mpd[!is.na(mpd)]
tu = data.frame(fam = tca$cat2[idxs], mpd = mpd, stringsAsFactors = F)
tu = tu[tu$fam %in% fams,]

breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
tu = cbind(tu, dis = cut(tu$mpd, breaks, include.lowest = T))
to = ddply(tu, .(fam, dis), summarise, cnt = length(fam))
to2 = ddply(tu, .(fam), summarise, cnt_fam = length(fam))
to = merge(to, to2, by = 'fam')
to = cbind(to, pct = to$cnt / to$cnt_fam)

tos = to[to$dis == '[0,0.01]',]
tos = tos[order(tos$pct),]
fams = tos$fam
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
  theme(legend.position = "top", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(2,1,0,5), "lines")) +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 8, colour = "grey", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0)) +
  annotation_custom(grob = textGrob(label = "Gene Family | #", just = c('right', 'bottom'), gp = gpar(fontsize = 8)), ymin = -0.05, ymax = -0.05, xmin = length(fams)+1, xmax = length(fams)+1)

for (i in 1:length(fams)) {
  p <- p + annotation_custom(grob = textGrob(label = fams[i], just = c('right'), gp = gpar(fontsize = 8)), ymin = -0.08, ymax = -0.08, xmin = i, xmax = i)
}
  
gt1 <- ggplot_gtable(ggplot_build(p))
gt1$layout$clip[gt1$layout$name == "panel"] <- "off"

fp = file.path(diro, "76.proteome.mpd.pdf")
pdf(file = fp, width = 6, height = 6, bg = 'transparent')
grid.newpage()
grid.draw(gt1)
dev.off()

##### gene fam AFS
fam_cnt = table(tca$cat2)
to0 = cbind(tca, n_org = norgs, stringsAsFactors = F)
to1 = ddply(to0, .(n_org, cat2), summarise, cnt = length(cat2))
to2 = cbind(to1, fam_cnt = fam_cnt[to1$cat2])
to = cbind(to2, fam_pct = to2$cnt / to2$fam_cnt)
to = to[to$cat2 %in% fams,]

tos = to[to$n_org == 16,]
fams1 = tos$cat2[order(tos$fam_pct)]
all_fams = unique(to$cat2)
fams1 = c(fams1, all_fams[!all_fams %in% fams1])
to$cat2 = factor(to$cat2, levels = fams1)
cnts = format(fam_cnt[fams1], big.mark = ',')

p2 = ggplot(to, aes(x = cat2, y = fam_pct, fill = n_org)) +
  geom_bar(stat = 'identity', position = "stack", width = 0.7) + 
#  scale_fill_manual(name = "AA distance:", breaks = labs, labels = labs, values = cols, guide = guide_legend(nrow = 2, byrow = T, label.position = "right", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  labs(fill = "Sharing Accessions") +
  coord_flip() +
  scale_x_discrete(name = '', labels = fams1) +
  scale_y_continuous(name = '', expand = c(0, 0)) +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.1, 'lines')) +
  theme(legend.position = "top", legend.key.size = unit(0.6, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,1,0,0), "lines")) +
  theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 8, colour = "black")) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = sprintf("%s/76.proteome.afs.pdf", diro)
ggsave(p2, filename = fp, width = 5, height = 5)


orgs_in = get_orgs( opt = "ingroup" )
nsams = apply(tst[,orgs_in], 1, function(x) sum(x != '3'))
ndels = apply(tst[,orgs_in], 1, function(x) sum(x %in% c('-', '')))
pis = 2*(ndels/nsams)*(1-ndels/nsams)

idxs = which(nsams > 5)
ti = cbind(tca[idxs,], ndel = ndels[idxs], stringsAsFactors = F)
tis = ti[ti$cat2 %in% fams,]
to = ddply(tis, .(cat2), summarise, cnt = length(cat2), prop = sum(ndel>0)/cnt)

ti = cbind(tca[idxs,], pi = pis[idxs], stringsAsFactors = F)
tis = ti[ti$cat2 %in% fams,]
to = ddply(tis, .(cat2), summarise, cnt = length(cat2), q25 = quantile(pi, 0.25), q50 = quantile(pi, 0.5), q75 = quantile(pi, 0.75))


#

tracks = list(gt1, gt2)
wds = c(1,4,4)

gt <- gtable(widths = unit(wds, "null"), height = unit(1, "null"))

gt <- gtable_add_grob(gt, gt1[, 2:3], 1, 1)
gt <- gtable_add_grob(gt, gt1[, 4], 1, 2)
gt <- gtable_add_grob(gt, gt2[, 4], 1, 3)



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
