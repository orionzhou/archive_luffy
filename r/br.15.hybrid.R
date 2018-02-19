require(grid)
require(tidyverse)
require(gtable)
require(RColorBrewer)
require(viridis)
require(Hmisc)
require(pheatmap)
source("common.R")

dirg = file.path(Sys.getenv("genome"), "Zmays_v4")
dirw = file.path(Sys.getenv("misc2"), "briggs")
diro = file.path(dirw, "41.qc")

fi = file.path(dirw, '32.bytissue.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

### expressed genes in each tissue ------------
tr2 = ti %>%
	group_by(Tissue, Genotype) %>%
	summarise(nexp = sum(FPM >= 0.1)) %>%
	mutate(flt = "FPM >= 0.1")

tr3 = ti %>%
	group_by(Tissue, Genotype) %>%
	summarise(nexp = sum(FPM >= 1)) %>%
	mutate(flt = "FPM >= 1")

to = rbind(tr2, tr3)

to$Genotype = factor(to$Genotype, levels = unique(ti$Genotype))
to$Tissue = factor(to$Tissue, levels = rev(unique(ti$Tissue)))
to$flt = factor(to$flt, levels = unique(to$flt))

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
ggsave(p1, filename = fo, width = 6, height = 8)

### visual
to = to[to$flt == 'FPM >= 1',]
to = to[order(to$Tissue, to$Genotype),]
tissues = unique(ti$Tissue)
tiss_map = 1:length(tissues)
names(tiss_map) = tissues
to = cbind(x = 1:nrow(to), to)
tox = to %>%
	group_by(Tissue) %>%
	summarise(xmax = max(x), x = mean(x))

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
theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 9, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 9, colour = "black", angle = 0, hjust = 1))
fo = file.path(diro, "21.expressed.genes2.pdf")
ggsave(p1, filename = fo, width = 6, height = 9)


### SPE genes
tissues = unique(ti$Tissue)
tissues = tissues[! tissues %in% c("kernel_14DAP", "endosperm_27DAP")]
to1 = filter(ti, Tissue %in% tissues & Genotype != "MxB")
to2 = spread(to1[,-5], Genotype, FPM)
to3 = mutate(to2, type1 = B73 >= 1 & Mo17 >= 1 & BxM >= 1,
		type2 = B73 >= 1 & Mo17 < 1 & BxM >= 1,
		type3 = B73 < 1 & Mo17 >= 1 & BxM >= 1,
	)

to4 = to3 %>% group_by(Tissue) %>%
	summarise(spe_b = sum(type2), spe_m = sum(type3))
to4 = gather(to4, type, ngene, -Tissue)
to4$Tissue = factor(to4$Tissue, levels = tissues)

p1 = ggplot(to4) +
  geom_bar(aes(x = Tissue, y = ngene, fill = type), stat="identity", position='dodge', width=0.6) +
  #scale_x_discrete(name = '# Shared Tissues') +
  scale_y_continuous(name = '# Genes') +
  scale_fill_brewer(palette = 'Set1', labels = c("B73-specific", "Mo17-specific")) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 30, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = file.path(diro, "25.spe.pdf")
ggsave(p1, filename = fp, width = 8, height = 4)


to4 = to3 %>% 
	group_by(gid) %>%
	summarise(ntis1 = sum(type1), ntis2 = sum(type2), ntis3 = sum(type3))
to5 = gather(to4, type, ntis, -gid)
to6 = to5[to5$ntis > 0,]
p1 = ggplot(to6) +
  geom_histogram(aes(x = ntis), stat = "count") +
  scale_x_discrete(name = '# Shared Tissues') +
  #scale_y_continuous(name = 'F1 B73 Proportion', expand = c(0.01, 0)) +
  scale_fill_manual(name = '', values = cols) +
  facet_wrap( ~ type, nrow = 1, scale = 'free') +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = file.path(diro, "26.spe.tissue.pdf")
ggsave(p1, filename = fp, width = 8, height = 3)


# heatmap
tk = to3[to3$type2 | to3$type3, c(1,2,7,8)]
tk = cbind(tk, type = 0)
tk$type[tk$type2] = 1
tk$type[tk$type3] = -1
tk = spread(tk[,c(1,2,5)], Tissue, type)
tk[is.na(tk)] = 0

hcl = hclust(dist(tk[,-1]), method = "ward.D")
gidsO = tk$gid[hcl$order]

tkl = gather(tk, tissue, type, -gid)
tkl$gid = factor(tkl$gid, levels = gidsO)
tkl$type = as.character(tkl$type)
tkl$tissue = factor(tkl$tissue, levels = tissues)

p1 = ggplot(tkl) +
  geom_tile(aes(x = tissue, y = gid, fill = type)) + 
  #scale_x_discrete(name = '') +
  scale_y_discrete(name = 'Genes') +
  scale_fill_manual(values = c("white", "#E41A1C", "#4DAF4A"), labels = c("non SPE", "SPE_B", "SPE_M")) + 
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill='grey', size=1)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 70, hjust = 1)) +
  theme(axis.text.y = element_blank())
fo = sprintf("%s/23.spe.heatmap.pdf", diro)
ggsave(p1, filename = fo, width = 4.5, height = 12)



