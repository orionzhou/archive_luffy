require(rtracklayer)
require(plyr)
require(grid)
require(ggplot2)

qnames1 = c(
  "HM058", "HM125", "HM056", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340"
)
dir = sprintf("%s/gene_sv", Sys.getenv("misc3"))

##### SV impact on gene families
qname = "HM004"
tname = "HM101"

tdir = sprintf("%s/%s", Sys.getenv("genome"), tname)
tfg = file.path(tdir, "51.tbl")
tg = read.table(tfg, sep = "\t", header = F, as.is = T)
colnames(tg) = c('chr', 'beg', 'end', 'srd', 'id', 'type', 'cat')
tg = tg[tg$type == 'cds',]
tg_size = ddply(tg, .(id, cat), summarise, size = sum(end - beg + 1))


qdir = sprintf("%s/%s", Sys.getenv("genome"), qname)
qfg = file.path(qdir, "51.tbl")
qg = read.table(qfg, sep = "\t", header = F, as.is = T)
colnames(qg) = c('chr', 'beg', 'end', 'srd', 'id', 'type', 'cat')
qg = qg[qg$type == 'cds',]
qg_size = ddply(qg, .(id, cat), summarise, size = sum(end - beg + 1))

do = data.frame()

type = 'ins'
f_bed = sprintf("%s/%s.ins.bed", dir, qname)
gene = qg
gene_size = qg_size
t_bed = read.table(f_bed, sep = "\t", header = F, as.is = T)[,c(1:3,7)]
colnames(t_bed) = c("chr", "beg", "end", "len")
t_bed$beg = t_bed$beg + 1
tm = merge(gene, t_bed, by = c("chr", "beg", "end"), all = T)
tms = ddply(tm, .(id), summarise, len = sum(len))
tp = merge(gene_size, tms, by = "id", all = T)
tp = cbind(tp, pct = tp$len / tp$size)
td = ddply(tp, .(cat), summarise, mean = mean(pct), std = sd(pct))
do = rbind(do, cbind(type = type, td))

type = 'gan'
f_bed = sprintf("%s/%s.gan.bed", dir, qname)
gene = qg
gene_size = qg_size
t_bed = read.table(f_bed, sep = "\t", header = F, as.is = T)[,c(1:3,7)]
colnames(t_bed) = c("chr", "beg", "end", "len")
t_bed$beg = t_bed$beg + 1
tm = merge(gene, t_bed, by = c("chr", "beg", "end"), all = T)
tms = ddply(tm, .(id), summarise, len = sum(len))
tp = merge(gene_size, tms, by = "id", all = T)
tp = cbind(tp, pct = tp$len / tp$size)
td = ddply(tp, .(cat), summarise, mean = mean(pct), std = sd(pct))
do = rbind(do, cbind(type = type, td))

type = 'del'
f_bed = sprintf("%s/%s.del.bed", dir, qname)
gene = tg
gene_size = tg_size
t_bed = read.table(f_bed, sep = "\t", header = F, as.is = T)[,c(1:3,7)]
colnames(t_bed) = c("chr", "beg", "end", "len")
t_bed$beg = t_bed$beg + 1
tm = merge(gene, t_bed, by = c("chr", "beg", "end"), all = T)
tms = ddply(tm, .(id), summarise, len = sum(len))
tp = merge(gene_size, tms, by = "id", all = T)
tp = cbind(tp, pct = tp$len / tp$size)
td = ddply(tp, .(cat), summarise, mean = mean(pct), std = sd(pct))
do = rbind(do, cbind(type = type, td))

type = 'los'
f_bed = sprintf("%s/%s.los.bed", dir, qname)
gene = tg
gene_size = tg_size
t_bed = read.table(f_bed, sep = "\t", header = F, as.is = T)[,c(1:3,7)]
colnames(t_bed) = c("chr", "beg", "end", "len")
t_bed$beg = t_bed$beg + 1
tm = merge(gene, t_bed, by = c("chr", "beg", "end"), all = T)
tms = ddply(tm, .(id), summarise, len = sum(len))
tp = merge(gene_size, tms, by = "id", all = T)
tp = cbind(tp, pct = tp$len / tp$size)
td = ddply(tp, .(cat), summarise, mean = mean(pct), std = sd(pct))
do = rbind(do, cbind(type = type, td))


f_tlg = sprintf("%s/%s.tlc.gan.bed", dir, qname)


f_tll = sprintf("%s/%s.tlc.los.bed", dir, qname)

fo = sprintf("%s/compstat/comp_genesv_%s.pdf", Sys.getenv("misc3"), qname)
p = ggplot(do) +
  geom_bar(mapping = aes(x = cat, y = mean), 
    stat = 'identity', geom_params=list(width = 0.8, alpha = 0.8)) +
  geom_errorbar(mapping = aes(x = cat, ymin = mean-std, ymax = mean+std), 
    stat = 'identity', geom_params=list(width = 0.4, alpha = 0.8)) +
#  scale_fill_brewer(palette='Set1', name = '', guide = guide_legend(nrow = 1, byrow = T, label.position = "right", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  coord_flip() +
  scale_x_discrete(name = 'fam') +
  scale_y_continuous(name = 'pct') +
  facet_wrap( ~ type, ncol = 4) +
  theme_bw() +
  theme(legend.position = "top", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(0,1,1,0), "lines")) +
  theme(axis.text.x = element_text(size = 7, colour = "gray", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0))
ggsave(p, filename = fo, width = 8, height = 10)


