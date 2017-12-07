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
diro = file.path(dirw, "44.hybrid")

fh = file.path(dirw, "00.1.read.correct.tsv")
th = read.table(fh, sep = "\t", header = T, as.is = T)[,1:5]
th$Genotype[th$Genotype == 'B73xMo17'] = 'BxM'
th$Genotype[th$Genotype == 'Mo17xB73'] = 'MxB'

### 
ff = file.path(dirw, '36.long.filtered.tsv')
tf = read.table(ff, header = T, sep = "\t", as.is = T)
tf = tf[tf$Genotype != 'Mo17xB73',]
tf$Genotype = factor(tf$Genotype, levels = c("B73","Mo17", "B73xMo17"))
tf$fpm = asinh(tf$fpm)

cv <- function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)*100
quantilecount <- function(x, fpm_min, fpm_max) table(cut(x, breaks = seq(from=fpm_min,to=fpm_max,length.out=5),include.lowest=T))

tf2 = spread(tf[,-5], Genotype, fpm)
grp = dplyr::group_by(tf2, gid)
tf3 = dplyr::summarise(grp, 
    fpm_min = min(B73, Mo17, B73xMo17),
    fpm_max = max(B73, Mo17, B73xMo17),
    cvb = sd(B73)/mean(B73)*100, 
    cvm = sd(Mo17)/mean(Mo17)*100,
    bq1 = quantilecount(B73, fpm_min, fpm_max)[1],
    bq23 = sum(quantilecount(B73, fpm_min, fpm_max)[2:3]),
    bq4 = quantilecount(B73, fpm_min, fpm_max)[4],
    mq1 = quantilecount(Mo17, fpm_min, fpm_max)[1],
    mq23 = sum(quantilecount(Mo17, fpm_min, fpm_max)[2:3]),
    mq4 = quantilecount(Mo17, fpm_min, fpm_max)[4],
    hq1 = quantilecount(B73xMo17, fpm_min, fpm_max)[1],
    hq23 = sum(quantilecount(B73xMo17, fpm_min, fpm_max)[2:3]),
    hq4 = quantilecount(B73xMo17, fpm_min, fpm_max)[4]
)
tf4 = filter(tf3, bq1 >=3 & bq4 >= 3 & bq23 >= 6 &
    mq1 >= 3 & mq4 >= 3 & mq23 >= 6 &
    hq1 >= 3 & hq4 >= 3 & hq23 <= 2)
nrow(tf4)


gid = tf4$gid[1]
for (gid in tf4$gid) {
tp = tf[tf$gid == gid,]
tpb = tp[tp$Genotype == "B73xMo17",]
tissues = tpb$Tissue[order(tpb$fpm)]
tp$Tissue = factor(tp$Tissue, levels = tissues)

ys = seq(from=min(tp$fpm),to=max(tp$fpm),length.out=5)

p1 = ggplot(tp) +
  geom_point(aes(x = Tissue, y = fpm, color = Genotype, shape = Genotype), size = 2) + 
  geom_line(aes(x = Tissue, y = fpm, color = Genotype, group = Genotype),) + 
  #scale_x_continuous(breaks = tox$x, labels = tox$Tissue, limits = c(0,max(to$x)+1), expand = c(0,0)) +
  scale_y_continuous(name = 'asinh(FPM)', breaks = ys, expand = c(0,0)) +
  #geom_hline(yintercept = ys, color='grey90') +
  scale_color_manual(values = brewer.pal(3,"Set1"), name = 'Genotype:') +
  scale_shape_manual(values = c(16,15,17), name = 'Genotype:') +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = c(0.1,0.9), legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_text(size=8), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 30, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fo = sprintf("%s/21.%s.pdf", diro, gid)
ggsave(p1, filename = fo, width = 6, height = 5)
}




tf3 = dplyr::summarise(grp, bmc = cor(B73,Mo17), 
    pcc.b = cor(B73[order(B73)], 1:17),
    pcc.m = cor(Mo17[order(Mo17)], 1:17),
    #pcc.h = cor(B73xMo17[order(B73)], 1:17),
    pval.b = cor.test(B73[order(B73)], 1:17)$p.value,
    pval.m = cor.test(Mo17[order(Mo17)], 1:17)$p.value,
    #pval.h = cor.test(B73xMo17[order(B73)], 1:17)$p.value
)


tf3 = tf[tf$gid %in% gids_bms,]
tf3 = spread(tf3[,-5], Genotype, fpm)

grp = dplyr::group_by(tf3, gid)
tf4 = dplyr::summarise(grp, 
)

tf5 = data.frame(tf4[tf4$pcc.b > 0.8 & tf4$pcc.m > 0.8,])