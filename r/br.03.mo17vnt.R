require(plyr)
require(dplyr)
require(ggplot2)

dirw = '/home/springer/zhoux379/scratch/mo17vnt'
fi = file.path(dirw, "52.vnt/Mo17.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
tf = ti[ti$PASS==1,]



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

fr = file.path(dirw, "stats/00.2.rpc.pdf")
ggsave(pr, filename=fr, width=5, height=6)
