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

# 
fi = file.path(dirw, "01.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)


# Reg based similarity
ti2 = ti[!ti$Tissue %in% c("endosperm_14DAP", "endosperm_27DAP", "kernel_14DAP"),]

grp = group_by(ti2, gid)
tig = summarise(grp, 
	ntis_a = sum(!is.na(Reg)), 
	ntis_reg = sum(is.de != 'non-DE' & !is.na(Reg)),
	ntis_cis = sum(!is.na(Reg) & Reg == 'cis'),
	ntis_trans = sum(!is.na(Reg) & Reg == 'trans')
)
gids = tig$gid[tig$ntis_a >= 10 & tig$ntis_reg >= 1]

ti3 = ti2[ti2$gid %in% gids, c("Tissue","gid","Reg")]
#ti3$Reg[!is.na(ti3$Reg) & ti3$Reg == 'cis+trans'] = 'none'
ti3$Reg[is.na(ti3$Reg)] = 'none'
ti3$Reg[!ti3$Reg %in% c("cis", "trans")] = 'none'

#reg_map = c("Mo17-biased"=-2, 'trans'=-1, 'none'=0, 'cis'=1, 'B73-biased'=2)
reg_map = c('trans'=-1, 'none'=0, 'cis'=1)
ti3 = cbind(ti3, Reg2 = reg_map[ti3$Reg])
ti4 = spread(ti3[,-3], Tissue, Reg2)

tk = ti4
#gower_dist = daisy(tk[,-1], metric = "gower")

hc = hclust(dist(tk[,-1]), method = "ward.D")
gidsO = tk$gid[hc$order]

tkl = ti3
tkl$gid = factor(tkl$gid, levels = gidsO)
tkl$Tissue = factor(tkl$Tissue, levels = unique(tl$Tissue))
tkl$Reg = factor(tkl$Reg, levels = names(reg_map))
cols = brewer.pal(4, "Set1")
cols = c(cols[1:2], "white", cols[3:4])
cols = c(cols[1], "white", cols[2])

p1 = ggplot(tkl) +
  geom_tile(aes(x = Tissue, y = gid, fill = Reg)) + 
  #scale_x_discrete(name = '') +
  scale_y_discrete(name = 'Genes') +
  scale_fill_manual(values = cols) + 
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill='grey', size=1)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 70, hjust = 1)) +
  theme(axis.text.y = element_blank())
fo = sprintf("%s/25.heatmap.reg.pdf", dirw)
ggsave(p1, filename = fo, width = 4, height = 8)

# sharing pattern of cis and trans across tissues
tig2 = tig[tig$gid %in% gids,]
tp1 = data.frame(type = 'cis', table(tig2$ntis_cis), stringsAsFactors = F)
tp2 = data.frame(type = 'trans', table(tig2$ntis_trans), stringsAsFactors = F)
tp = rbind(tp1, tp2)
colnames(tp) = c("type", "nTissue", "nGene")
tp = tp[tp$nTissue != 0,]

p1 = ggplot(tp) +
  geom_histogram(aes(x = nTissue, y = nGene, fill = type), position = 'dodge', width = 0.7, stat = 'identity') +  	
  scale_x_discrete(name = '# Sharing Tissues') +
  scale_y_continuous(name = '# Genes') +
  scale_fill_brewer(palette = "Set1") + 
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = c(0.8, 0.9), legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fo = sprintf("%s/30.sharing.reg.pdf", dirw)
ggsave(p1, filename = fo, width = 4, height = 4)

## non-DE DOA/HoM
ti2 = ti[!ti$Tissue %in% c("endosperm_14DAP", "endosperm_27DAP", "kernel_14DAP"),]

ti3 = ti2[ti2$is.de == 'non-DE' & !is.na(ti2$HoM) & (ti2$B73 >=1 | ti2$Mo17 >=1),]

describe(ti3$HoM)
ti4 = ti3
ti4$HoM[!is.na(ti4$HoM) & ti4$HoM < -2] = -2
ti4$HoM[!is.na(ti4$HoM) & ti4$HoM > 2] = 2

p1 = ggplot(ti4) +
  geom_density(aes(x = HoM)) + 
  scale_x_continuous(name = 'log(F1/MP) (non-DE genes)') +
  #scale_fill_brewer(palette = "Set1") + 
  facet_wrap(~Tissue, nrow=4) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8))
fo = sprintf("%s/19.HoM.pdf", dirw)
ggsave(p1, filename = fo, width = 6, height = 6)
