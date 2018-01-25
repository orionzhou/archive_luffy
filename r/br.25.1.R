require(plyr)
require(dplyr)
require(tidyr)
require(ape)
require(ggplot2)
require(Hmisc)
require(RColorBrewer)

dirw = file.path(Sys.getenv("misc2"), "briggs", "44.hybrid")

fl = file.path(dirw, '../00.1.read.correct.tsv')
tl = read.table(fl, header = T, sep = "\t", as.is = T)[,1:5]
tl = within(tl, {label = sprintf("%s: %s %s rep%d", SampleID, Tissue, Genotype, Treatment)})

### build master table
fi = file.path(dirw, "../35.long.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

ti2 = spread(ti[,-4], Genotype, fpkm)
ti3 = within(ti2, {
    BoM = log(B73 / Mo17)
    MP = (B73 + Mo17) / 2
    DoA = (B73xMo17 - MP) / (pmax(B73, Mo17) - MP)
    HoM = log(B73xMo17 / MP)
    rm(Mo17xB73, MP)
})
ti4 = ti3[,c(1:5,8:6)]

head(ti4)
nrow(ti4)/39005

describe(ti4$BoM)
ti4$BoM[is.nan(ti4$BoM)] = NA
ti4$BoM[ti4$BoM == -Inf] = -9
ti4$BoM[ti4$BoM == Inf] = 9
describe(ti4$BoM)

describe(ti4$DoA)
ti4$DoA[is.nan(ti4$DoA)] = NA
ti4$DoA[ti4$DoA == Inf] = 70000
describe(ti4$DoA)

describe(ti4$HoM)
ti4$HoM[is.nan(ti4$HoM)] = NA
ti4$HoM[ti4$HoM == -Inf] = -7
ti4$HoM[ti4$HoM == Inf] = 7
describe(ti4$HoM)

# add DE
ff = file.path(dirw, '../43.deseq/01.de.tsv')
tf = read.table(ff, header = T, sep = "\t", as.is = T)
tf = tf[tf$comp == 'B73 vs Mo17', -2]
colnames(tf)[1] = "Tissue"

ti5 = merge(ti4, tf, by = c("Tissue", "gid"), all.x = T)
ti5$is.de[is.na(ti5$is.de)] = 'non-DE'
nrow(ti5)/39005

# add ASE
fb = file.path(dirw, "../46.ase/30.ase.tsv")
tb = read.table(fb, sep = "\t", header = T, as.is = T)

ti6 = merge(ti5, tb, by = c("Tissue", "gid"), all.x = T)
nrow(ti6)/39005

to = ti6
fo = file.path(dirw, "01.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)


### 
fi = file.path(dirw, "01.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)


# FPKM based tissue similarity
grp = group_by(ti, gid)
tog = summarise(grp, ntis_b = sum(B73>=1), ntis_m = sum(Mo17>=1))
gids = tog$gid[tog$ntis_b > 0 & tog$ntis_m > 0]

to = ti
to2 = to[to$gid %in% gids, c("Tissue","gid","B73","Mo17")]
to3 = gather(to2, gt, FPKM, -Tissue, -gid)
to3 = within(to3, {
	tis_gt = sprintf("%s|%s", Tissue, gt)
	rm(Tissue, gt)
})
to4 = spread(to3, tis_gt, FPKM)

e = asinh(to4[,-1])
pca <- prcomp(e, center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
tp = cbind.data.frame(tis_gt = rownames(x), x[,1:5], stringsAsFactors = F)
res = strsplit(tp$tis_gt, split = "[|]")
tp = cbind(tp, Tissue = sapply(res, "[", 1), Genotype = sapply(res, "[", 2))
tp$Tissue = factor(tp$Tissue, levels = unique(tl$Tissue))
tp$Genotype = factor(tp$Genotype, levels = unique(tp$Genotype))

cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

p1 = ggplot(tp) +
  geom_point(aes(x = PC1, y = PC2, shape = Genotype, color = Tissue), size = 4) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(name = ylab) +
  scale_color_manual(name = "", values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.3,0.1,0.1,0.1), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fp = sprintf("%s/10.pca.FPKM.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 5.5)


# BoM based similarity
to = ti
grp = group_by(to, gid)
tog = summarise(grp, ntis_de = sum(is.de != 'non-DE'))
gids = tog$gid[tog$ntis_de >= 1]

to2 = to[to$gid %in% gids, c("Tissue","gid","BoM")]
to3 = spread(to2, Tissue, BoM)
to3[is.na(to3)] = 0

e = to3[,-1]
pca <- prcomp(e, center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
tp = cbind.data.frame(Tissue = rownames(x), x[,1:5], stringsAsFactors = F)
tp$Tissue = factor(tp$Tissue, levels = unique(tl$Tissue))

cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

p1 = ggplot(tp) +
  geom_point(aes(x = PC1, y = PC2, color = Tissue), size = 4) +
  geom_text(aes(x = PC1, y = PC2, label = Tissue), size = 3, nudge_y = 0.03, check_overlap = T) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(name = ylab) +
  scale_color_manual(name = "", values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fp = sprintf("%s/11.pca.BoM.pdf", dirw)
ggsave(p1, filename = fp, width = 5.5, height = 5.5)

# BoM heatmap
tk = to3
hc = hclust(dist(tk[,-1]), method = "ward.D")
gidsO = tk$gid[hc$order]

tkl = gather(tk, tissue, BoM, -gid)
tkl$gid = factor(tkl$gid, levels = gidsO)
tkl$tissue = factor(tkl$tissue, levels = unique(tl$Tissue))

p1 = ggplot(tkl) +
  geom_tile(aes(x = tissue, y = gid, fill = BoM)) + 
  #scale_x_discrete(name = '') +
  scale_y_discrete(name = 'Genes') +
  scale_fill_gradient2() + #(colours = terrain.colors(10)) + 
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill='grey', size=1)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 70, hjust = 1)) +
  theme(axis.text.y = element_blank())
fo = sprintf("%s/21.heatmap.BoM.pdf", dirw)
ggsave(p1, filename = fo, width = 4.5, height = 12)


## DoA similarity
ti2 = ti[!ti$Tissue %in% c("endosperm_14DAP", "endosperm_27DAP", "kernel_14DAP"),]

grp = group_by(ti2, gid)
tig = summarise(grp, 
	ntis_de = sum(is.de != 'non-DE')
)
gids = tig$gid[tig$ntis_de >= 10]

ti2 = ti2[ti2$gid %in% gids, c("Tissue","gid","is.de","DoA")]
ti2$DoA[ti2$is.de == 'non-DE'] = NA

describe(ti2$DoA)
ti2$DoA[!is.na(ti2$DoA) & ti2$DoA < -3] = -3
ti2$DoA[!is.na(ti2$DoA) & ti2$DoA > 3] = 3

ti3 = spread(ti2[,-3], Tissue, DoA)

e = ti3[,-1]
e[is.na(e)] = 0

pca <- prcomp(e, center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
tp = cbind.data.frame(Tissue = rownames(x), x[,1:5], stringsAsFactors = F)
tp$Tissue = factor(tp$Tissue, levels = unique(tl$Tissue))

cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

p1 = ggplot(tp) +
  geom_point(aes(x = PC1, y = PC2, color = Tissue), size = 4) +
  geom_text(aes(x = PC1, y = PC2, label = Tissue), size = 3, nudge_y = 0.03, check_overlap = T) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(name = ylab) +
  scale_color_manual(name = "", values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fp = sprintf("%s/12.pca.DoA.pdf", dirw)
ggsave(p1, filename = fp, width = 5.5, height = 5.5)


# DoA heatmap
tk = ti3
tk[is.na(tk)] = 0
hc = hclust(dist(tk[,-1]), method = "ward.D")
gidsO = tk$gid[hc$order]

tkl = ti2
tkl$gid = factor(tkl$gid, levels = gidsO)
tkl$Tissue = factor(tkl$Tissue, levels = unique(tl$Tissue))

p1 = ggplot(tkl) +
  geom_tile(aes(x = Tissue, y = gid, fill = DoA)) + 
  #scale_x_discrete(name = '') +
  scale_y_discrete(name = 'Genes') +
  scale_fill_gradient2(breaks=c(-2,0,2), labels=c("F1<MP", "F1=MP", "F1>HP")) + #(colours = terrain.colors(10)) + 
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill='grey', size=1)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 70, hjust = 1)) +
  theme(axis.text.y = element_blank())
fo = sprintf("%s/22.heatmap.DoA.pdf", dirw)
ggsave(p1, filename = fo, width = 4, height = 10)

# DOA - DE
ti2 = ti[!ti$Tissue %in% c("endosperm_14DAP", "endosperm_27DAP", "kernel_14DAP"),]

grp = group_by(ti2, gid)
tig = summarise(grp, 
	ntis_de = sum(is.de != 'non-DE')
)
gids = tig$gid[tig$ntis_de >= 1]

ti3 = ti2[ti2$gid %in% gids, c("Tissue","gid","p.propB","is.de","DoA")]
ti4 = ti3[ti3$is.de != 'non-DE' & !is.na(ti3$p.propB),]

describe(ti4$DoA)
ti4$DoA[!is.na(ti4$DoA) & ti4$DoA < -2] = -2
ti4$DoA[!is.na(ti4$DoA) & ti4$DoA > 2] = 2

p1 = ggplot(ti4) +
  geom_point(aes(x = DoA, y = p.propB), size = 0.5) + 
  scale_x_continuous(name = 'DoA: (F1-MP)/(HP-MP)') +
  scale_y_continuous(name = 'Parental B73 Proportion') +
  #scale_fill_brewer(palette = "Set1") + 
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8))
fo = sprintf("%s/51.DoA.DE.pdf", dirw)
ggsave(p1, filename = fo, width = 5, height = 5)

# DoA - AR
ti2 = ti[!ti$Tissue %in% c("endosperm_14DAP", "endosperm_27DAP", "kernel_14DAP"),]

grp = group_by(ti2, gid)
tig = summarise(grp, 
	ntis_de = sum(is.de != 'non-DE')
)
gids = tig$gid[tig$ntis_de >= 1]

ti3 = ti2[ti2$gid %in% gids, c("Tissue","gid","h.propB","is.de","DoA")]
ti4 = ti3[ti3$is.de != 'non-DE' & !is.na(ti3$h.propB),]

describe(ti4$DoA)
ti4$DoA[!is.na(ti4$DoA) & ti4$DoA < -2] = -2
ti4$DoA[!is.na(ti4$DoA) & ti4$DoA > 2] = 2

p1 = ggplot(ti4) +
  geom_point(aes(x = DoA, y = h.propB), size = 0.5) + 
  scale_x_continuous(name = 'DoA: (F1-MP)/(HP-MP)') +
  scale_y_continuous(name = 'Hybrid B73 Proportion') +
  #scale_fill_brewer(palette = "Set1") + 
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8))
fo = sprintf("%s/52.DE.DoA.pdf", dirw)
ggsave(p1, filename = fo, width = 5, height = 5)


## AR similarity
ti2 = ti[!is.na(ti$h.propB), c("Tissue","gid","h.propB")]
ti2 = ti2[!ti2$Tissue %in% c("endosperm_14DAP", "endosperm_27DAP", "kernel_14DAP"),]

ti3 = spread(ti2, Tissue, h.propB)

e = ti3[,-1]
e[is.na(e)] = 0.5
pca <- prcomp(e, center = F, scale. = F, na.action = 'na.omit')
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
tp = cbind.data.frame(Tissue = rownames(x), x[,1:5], stringsAsFactors = F)
tp$Tissue = factor(tp$Tissue, levels = unique(tl$Tissue))

cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

p1 = ggplot(tp) +
  geom_point(aes(x = PC1, y = PC2, color = Tissue), size = 4) +
  geom_text(aes(x = PC1, y = PC2, label = Tissue), size = 3, nudge_y = 0.03, check_overlap = T) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(name = ylab) +
  scale_color_manual(name = "", values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fp = sprintf("%s/13.pca.AR.pdf", dirw)
ggsave(p1, filename = fp, width = 5.5, height = 5.5)


# AR heatmap
tk = ti3
hc = hclust(dist(tk[,-1]), method = "ward.D")
gidsO = tk$gid[hc$order]

tkl = gather(tk, tissue, BoM, -gid)
tkl$gid = factor(tkl$gid, levels = gidsO)
tkl$tissue = factor(tkl$tissue, levels = unique(tl$Tissue))

p1 = ggplot(tkl) +
  geom_tile(aes(x = tissue, y = gid, fill = BoM)) + 
  #scale_x_discrete(name = '') +
  scale_y_discrete(name = 'Genes') +
  scale_fill_gradient2(name = '', midpoint = 0.5, na.value = 'grey50') +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill='grey', size=1)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 70, hjust = 1)) +
  theme(axis.text.y = element_blank())
fo = sprintf("%s/23.heatmap.AR.pdf", dirw)
ggsave(p1, filename = fo, width = 4.5, height = 12)

# DE - AR
ti2 = ti[!ti$Tissue %in% c("endosperm_14DAP", "endosperm_27DAP", "kernel_14DAP"),]
ti3 = ti2[!is.na(ti2$h.propB), c("Tissue","gid","BoM","is.de","p.propB","h.propB")]

p1 = ggplot(ti3) +
  geom_point(aes(x = h.propB, y = p.propB), size = 0.5) + 
  scale_x_continuous(name = 'F1 B73 Proportion') +
  scale_y_continuous(name = 'Parental B73 Proportion') +
  #scale_fill_brewer(palette = "Set1") + 
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8))
fo = sprintf("%s/52.DE.AR.pdf", dirw)
ggsave(p1, filename = fo, width = 5, height = 5)



