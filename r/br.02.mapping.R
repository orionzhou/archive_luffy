require(plyr)
require(dplyr)
require(tidyr)
require(ggplot2)

dirw = '/home/springer/zhoux379/scratch/briggs2'
dirw = file.path(Sys.getenv("misc2"), "briggs")
diro = file.path(dirw, "41.qc")

## raw read counts
fi = file.path(dirw, '00.2.trim.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)[,c(1:5,8,13)]

levs_tissue = unique(ti$Tissue)
levs_genotype = unique(ti$Genotype)
ti$Tissue = factor(ti$Tissue, levels = levs_tissue)
ti$Genotype = factor(ti$Genotype, levels = levs_genotype)
ti = ti[order(ti$Tissue, ti$Genotype),]

ti2 = cbind(x = 1:nrow(ti), ti[,c(2:7)])
tix = ddply(ti2, .(Tissue, Genotype), summarise, x = (min(x)+max(x))/2, lab = sprintf("%s_%s", Genotype[1], Tissue[1]))


pr = ggplot(ti2) +
  geom_rect(aes(xmin=x-0.4, xmax=x+0.4, ymin=0, ymax=RetainedReadPairCount/1000000, fill=Genotype), stat='identity') +
  geom_rect(aes(xmin=x-0.4, xmax=x+0.4, ymin=RetainedReadPairCount/1000000, ymax=ReadPairCount/1000000), fill='Black', stat='identity') +
#  layer(data=ti2, geom="bar", mapping=aes(x=Tissue, y=total/1000000 - 5, label=sprintf("%.01fx", cov_mapping)), geom_params=list(size=3.5, color='white')) +
  scale_fill_brewer(palette='Dark2') +
  scale_x_continuous(name='Tissue', breaks = tix$x, labels = tix$Tissue, expand = c(0,0), limits = c(0,max(ti2$x)+1)) +
  scale_y_continuous(name='Million Read Pairs (trimmed)', expand = c(0.01,0)) +
  coord_flip() +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.3,1), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1)) + 
  geom_rect(xmin=64,xmax=65,ymin=40,ymax=41,fill='black') +
  geom_label(x=64, y=41.5, hjust=0, vjust=0, label='Filtered Reads', label.padding = unit(0.2, "lines"), label.r = unit(0, "lines"), label.size = 0, size=2.5)

fr = file.path(diro, "00.2.rpc.pdf")
ggsave(pr, filename=fr, width=6, height=8)

### mapping stats
fi = file.path(dirw, "00.3.hisat.tsv")
ti = read.table(fi, header = T, sep = "\t", stringsAsFactors = F)[,c(1,3,4,17:24)]
levs_tissue = unique(ti$Tissue)
levs_genotype = unique(ti$Genotype)
ti$Tissue = factor(ti$Tissue, levels = levs_tissue)
ti$Genotype = factor(ti$Genotype, levels = levs_genotype)
ti = ti[order(ti$Tissue, ti$Genotype),]

ti2 = cbind(x = 1:nrow(ti), ti[,c(2:3,6:8)], Pair_Map_Lq = ti$Pair_Map-ti$Pair_Map_Hq)
to = gather(ti2, key = type, value = ReadPair, -x, -Tissue, -Genotype)
levs_type = c("Pair_Unmap", "Pair_Orphan", "Pair_Map_Lq", "Pair_Map_Hq")
to$type = factor(to$type, levels = levs_type)

tiss_gts = unique(sprintf("%s_%s", to$Tissue, to$Genotype))
tmap = 1:length(tiss_gts)
names(tmap) = tiss_gts
to = cbind(to, x2 = tmap[sprintf("%s_%s", to$Tissue, to$Genotype)])
to$x = to$x + (to$x2 - 1)

tox = ddply(to, .(Tissue, Genotype), summarise, x = (min(x)+max(x))/2, lab = sprintf("%s_%s", Genotype[1], Tissue[1]))

pr = ggplot(to) +
  geom_bar(aes(x=x, y=ReadPair/1000000, fill=type), stat='identity', position='stack') +
  scale_fill_brewer(palette='Set2') +
  scale_x_continuous(name='Tissue', breaks = tox$x, labels = tox$lab, expand = c(0,0)) +
  scale_y_continuous(name='Million Read Pairs', expand = c(0,0)) +
  coord_flip() +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.5,0.1,0.1), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.3,1), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))

fr = file.path(diro, "00.4.mapping.pdf")
ggsave(pr, filename=fr, width=6, height=10)

### benchmarking hisat2
fh = file.path(dirw, "../bt/00.3.hisat.tsv")
th = read.table(fh, header = T, sep = "\t", stringsAsFactors = F)[,c(1,3,4,17:24)]

tip = cbind(ti[,c(1:3)], pctmapv = (ti$Pair_Map+ti$Pair_Orphan/2)/ti$Pair*100)
thp = cbind(th[,c(1:3)], pctmapr = (th$Pair_Map+th$Pair_Orphan/2)/th$Pair*100)
tp = merge(thp[,c(1,4)], tip, by = 'SampleID')
tp$Tissue = factor(tp$Tissue, levels = unique(tp$Tissue))
tp$Genotype = factor(tp$Genotype, levels = unique(tp$Genotype))
tp = tp[order(tp$Tissue, tp$Genotype), c(3,4,2,5)]
tp = cbind(x = 1:nrow(tp), tp)

to = reshape(tp, direction='long', varying=list(4:5), idvar='x', timevar='type', v.names='pctmap', times = c("Hisat2 on B73", "Hisat2 on B73 + Mo17-variants"))
tox = ddply(to, .(Tissue, Genotype), summarise, minx = min(x), maxx = max(x), x = (min(x)+max(x))/2, lab = sprintf("%s_%s", Genotype[1], Tissue[1]))

pr = ggplot(to) +
  geom_bar(aes(x=x, y=pctmap, fill=type), stat='identity', position='dodge', width=0.7) +
  layer(data=tox, geom="segment", mapping=aes(x=minx, xend=maxx, y=0, yend=0), stat = 'identity', position = 'stack') +
  scale_fill_brewer(palette='Pastel1') +
  scale_x_continuous(name='Tissue', breaks = tox$x, labels = tox$lab, expand = c(0,0)) +
  scale_y_continuous(name='Reads Mapping Rate', expand = c(0.03,0), limits = c(0,100)) +
  coord_flip() +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.5,0.1,0.1), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.3,1), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))

fr = file.path(diro, "00.5.benchmark.hisat.pdf")
ggsave(pr, filename=fr, width=5, height=7)
