require(plyr)
require(ape)
require(tidyr)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(GenomicRanges)
require(pheatmap)
source("common.R")

dirg = file.path(Sys.getenv("genome"), "Zmays_v4")
dirw = file.path(Sys.getenv("misc2"), "briggs")
diro = file.path(dirw, "41.qc")

fh = file.path(dirw, "00.1.read.correct.tsv")
th = read.table(fh, sep = "\t", header = T, as.is = T)[,1:5]
th$Genotype[th$Genotype == 'B73xMo17'] = 'BxM'
th$Genotype[th$Genotype == 'Mo17xB73'] = 'MxB'

### 
fr = file.path(dirw, "32.rc.tsv")
tr = read.table(fr, header = T, sep = "\t", as.is = T)
 
ff = file.path(dirw, '33.fpm.tsv')
tf = read.table(ff, header = T, sep = "\t", as.is = T)

fk = file.path(dirw, '34.fpkm.tsv')
tk = read.table(fk, header = T, sep = "\t", as.is = T)

### expressed genes in each tissue
tr2 = gather(tr, sid, rc, -gid)
tr3 = merge(tr2, tm[,-2], by.x = 'sid', by.y = 'SampleID')

grp = dplyr::group_by(tr3, Genotype, Tissue, gid)
tr4 = dplyr::summarise(grp, rc = sum(rc))

grp = dplyr::group_by(tr4, Genotype, Tissue)
tt = dplyr::summarise(tr4, trc = sum(rc))
tr4 = merge(tr4, tt, by = c("Genotype", "Tissue"))
tr4 = cbind(tr4, fpm = tr4$rc/tr4$trc * 1000000)

to = data.frame()

grp = dplyr::group_by(tr4, Genotype, Tissue)

tr5 = dplyr::summarise(grp, nexp = sum(rc >= 1))
tr5 = cbind(flt = 'Read Count >= 1', data.frame(tr5))
to = rbind(to, tr5)

tr5 = dplyr::summarise(grp, nexp = sum(rc >= 2))
tr5 = cbind(flt = 'Read Count >= 2', data.frame(tr5))
to = rbind(to, tr5)

tr5 = dplyr::summarise(grp, nexp = sum(fpm >= 0.01))
tr5 = cbind(flt = 'FPM >= 0.01', data.frame(tr5))
to = rbind(to, tr5)

tr5 = dplyr::summarise(grp, nexp = sum(fpm >= 0.1))
tr5 = cbind(flt = 'FPM >= 0.1', data.frame(tr5))
to = rbind(to, tr5)

to$Genotype = factor(to$Genotype, levels = unique(tm$Genotype))
to$Tissue = factor(to$Tissue, levels = rev(unique(tm$Tissue)))
to$flt = factor(to$flt, levels = unique(to$flt))

#tr6 = spread(tr5, Genotype, nexp)
#tr6 = data.frame(tr6)

p1 = ggplot(to) +
  geom_tile(aes(x = Genotype, y = Tissue, fill = nexp)) + 
  geom_text(aes(x = Genotype, y = Tissue, label = nexp)) + 
  #scale_x_discrete(name = '') +
  #scale_y_discrete(name = 'RPKM') +
  scale_fill_gradientn(colours = brewer.pal(9, "Greens")) + 
  facet_wrap(~flt, nrow = 1) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype = 0)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'None') +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fo = file.path(diro, "21.expressed.genes.pdf")
ggsave(p1, filename = fo, width = 9, height = 8)


### visual
to = to[to$flt == 'FPM >= 0.1',]
to = to[order(to$Tissue, to$Genotype),]
tissues = unique(tm$Tissue)
tiss_map = 1:length(tissues)
names(tiss_map) = tissues
to = cbind(x = 1:nrow(to), to)
tox = ddply(to, .(Tissue), summarise, xmax = max(x), x = mean(x))

p1 = ggplot(to) +
  geom_point(aes(x = x, y = nexp, color = Genotype, shape = Genotype), size=2) + 
  geom_vline(xintercept = tox$xmax[-nrow(tox)]+0.5, alpha=0.3) +
  scale_x_continuous(breaks = tox$x, labels = tox$Tissue, limits = c(0,max(to$x)+1), expand = c(0,0)) +
  scale_y_continuous(name = 'Number Expressed Genes (FPM >= 0.1)') +
  scale_color_manual(values = rep(brewer.pal(3,"Set1")[1:2],2), name = 'Genotype:') +
  scale_shape_manual(values = c(16,15,1,0), name = 'Genotype:') +
  coord_flip() +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.5,0.5), legend.title = element_text(size=8), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 9, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 9, colour = "black", angle = 0, hjust = 1))
fo = file.path(diro, "21.expressed.genes2.pdf")
ggsave(p1, filename = fo, width = 6, height = 9)

###
fi = file.path(dirw, '35.long.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

ti = ti[ti$Genotype != 'Mo17xB73',]
ti$Genotype[ti$Genotype == 'B73xMo17'] = 'BxM'

td = spread(ti[,-4], Genotype, fpkm)
#td = td[td$B73 > 0 | td$Mo17 > 0,]
td = within(td, {
	ptag = ifelse(B73 > Mo17, "B > M", "M > B")
	doa = (BxM - (B73+Mo17)/2) / (pmax(B73,Mo17)-(B73+Mo17)/2)
})

p1 = ggplot(td) +
  geom_density(aes(x = doa, fill = ptag), alpha = 0.5) + 
  geom_vline(xintercept = 0, color = 'forestgreen') +
  scale_x_continuous(name = 'DOA: (F1-MP)/(HP-MP)', limits=c(-8,8)) +
  #scale_y_discrete(name = 'RPKM') +
  scale_fill_brewer(name = 'Direction', palette = "Set1") + 
  facet_wrap(~Tissue, nrow = 5) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  #theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype = 0)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = c(0.6,0.06), legend.direction = "vertical", legend.justification = c(0,0), legend.title = element_text(size=9), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fo = file.path(diro, "22.doa.pdf")
ggsave(p1, filename = fo, width = 8, height = 10)


ff = file.path(dirw, '43.deseq/01.de.tsv')
tf = read.table(ff, header = T, sep = "\t", as.is = T)
tf = tf[tf$comp == 'B73 vs Mo17', -2]

tp = merge(td, tf, by.x = c("Tissue", "gid"), by.y = c("tissue", "gid"))
tp$Tissue = factor(tp$Tissue, levels = unique(th$Tissue))

cutoff_fc = 4
tp2 = tp[tp$B73/tp$Mo17 >= cutoff_fc | tp$Mo17/tp$B73 >= cutoff_fc,]
p1 = ggplot(tp2) +
  geom_density(aes(x = doa, fill = ptag), alpha = 0.5) + 
  geom_vline(xintercept = 0, color = 'forestgreen') +
  scale_x_continuous(name = 'DOA: (F1-MP)/(HP-MP)', limits=c(-4,4)) +
  #scale_y_discrete(name = 'RPKM') +
  scale_fill_brewer(name = 'Direction', palette = "Set1") + 
  facet_wrap(~Tissue, nrow = 5) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  #theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype = 0)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = c(0.6,0.06), legend.direction = "vertical", legend.justification = c(0,0), legend.title = element_text(size=9), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fo = sprintf("%s/22.doa.de.%d.pdf", diro, cutoff_fc)
ggsave(p1, filename = fo, width = 8, height = 10)

source("go.R")
tx1 = go_enrich(tp$gid[tp$Tissue=='internode_v12' & tp$doa < -0.8 & tp$is.de == 'up'])
tx2 = go_enrich(tp$gid[tp$Tissue=='tasselstem_0DAP' & tp$doa < -0.8 & tp$is.de == 'up'])
tx3 = go_enrich(tp$gid[tp$Tissue=='tassel_v12' & tp$doa < -0.8 & tp$is.de == 'down'])
tx4 = go_enrich(tp$gid[tp$Tissue=='flagleaf_0DAP' & tp$doa < -0.8 & tp$is.de == 'down'])


### DOA heatmap
tds = merge(td, tf, by.x = c("gid", "Tissue"), by.y = c("gid", "tissue"))
colnames(tds)[2] = 'tissue'
sum(is.nan(tds$doa))

#use gidsO from DE

tds$gid = factor(tds$gid, levels = gidsO)
tds$tissue = factor(tds$tissue, levels = unique(ti$Tissue))
tds$doa[tds$doa > 5] = 5

p1 = ggplot(tds) +
  geom_tile(aes(x = tissue, y = gid, fill = doa)) + 
  #scale_x_discrete(name = '') +
  scale_y_discrete(name = 'Genes') +
  scale_fill_gradient2() + #(colours = terrain.colors(10)) + 
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.key.width = unit(0.8, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 70, hjust = 1)) +
  theme(axis.text.y = element_blank())
fo = sprintf("%s/23.doa.heatmap.pdf", diro)
ggsave(p1, filename = fo, width = 4, height = 12)


