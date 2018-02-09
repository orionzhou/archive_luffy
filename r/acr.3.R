require(plyr)
require(ggplot2)
require(dplyr)
require(tidyr)
require(Hmisc)

dirw = '/home/springer/zhoux379/data/misc1/maize.acr'

### log2B/P
fi = file.path(dirw, '10.rpm.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

ti = ti[ti$genotype %in% c("B73", "PH207"),]
tr = ti %>%
    group_by(gid, genotype, tissue) %>% 
    summarise(rpm = mean(rpm))
tr2 = spread(tr, genotype, rpm)
tr2 = within(tr2, {
    log2BP = log2(B73/PH207)
})
describe(tr2$log2BP[tr2$B73>0 & tr2$PH207>0])
tr2$log2BP[tr2$B73 == 0 | (is.finite(tr2$log2BP) & tr2$log2BP < -10)] = -10
tr2$log2BP[tr2$PH207 == 0 | (is.finite(tr2$log2BP) & tr2$log2BP > 10)] = 10
describe(tr2$log2BP)

tr = tr2[,c(1,2,5)]

# read mapping
fi = file.path(dirw, "09.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)


ti2 = merge(ti, tr, by = 'gid')
ti3 = ti2[ti2$varType %in% c("Identical", "qIns100", "tIns100"),]

p1 = ggplot(ti3) +
  #geom_histogram(aes(x = bpr, fill = varType), position='dodge', alpha=0.5) + 
  geom_density(aes(x = log2BP, fill = varType), alpha=0.4, linetype='blank') + 
  #geom_vline(xintercept = 0, color = 'forestgreen', size=0.1) +
  scale_x_continuous(name = 'log2(B/P)', limits=c(-10,10)) +
  #scale_y_discrete(name = 'RPKM') +
  scale_fill_brewer(name = 'Direction', palette = "Set1") + 
  facet_wrap(~tissue, ncol=1) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  #theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype = 0)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fo = sprintf("%s/31.log2BP.pdf", dirw)
ggsave(p1, filename = fo, width = 6, height = 8)


