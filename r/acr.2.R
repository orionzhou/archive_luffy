require(plyr)
require(ggplot2)
require(dplyr)
require(tidyr)
require(Hmisc)

dirw = '/home/springer/zhoux379/data/misc1/jackie.acr'

### read mapping table
fm = file.path(dirw, "09.tsv")
tm = read.table(fm, sep = "\t", header = T, stringsAsFactors = F)

### read expression table
fi = '/home/springer/nosha003/wgbs_schmitz/ACR/ACR_expression.txt'
ti = read.table(fi, sep = "\t", header = T, stringsAsFactors = F)
ti2 = ti[,33:37]
colnames(ti2)[1] = 'gene'
#ti2[ti2==1e-6] = 0
ti2[ti2 < 0.1] = 0
ti2 = within(ti2, {
	bpr.root = log(B73_root_rpm/PH207_root_rpm)
	bpr.leaf = log(B73_leaf_rpm/PH207_leaf_rpm)
})
# look at range of log ratio: [-8,8]
describe(ti2$bpr.root)
describe(ti2$bpr.leaf)
ti2$bpr.root[ti2$bpr.root == Inf] = 10
ti2$bpr.root[ti2$bpr.root == -Inf] = -10
ti2$bpr.root[is.na(ti2$bpr.root)] = NA
ti2$bpr.leaf[ti2$bpr.leaf == Inf] = 10
ti2$bpr.leaf[ti2$bpr.leaf == -Inf] = -10
ti2$bpr.leaf[is.na(ti2$bpr.leaf)] = NA
describe(ti2$bpr.root)
describe(ti2$bpr.leaf)

# only take genes in synteny with PH207 and maps uniquely
tm2 = tm[tm$syntag == 'inSynteny' & tm$maptag == 'Unique',]
tm3 = merge(tm2, ti[,33:37], by.x = 'gene', by.y = 'transcript')

# group genes into 'Identical' where all ACRs associated with that gene is 'Identical'
# and 'Polymorphic' if any of the associated ACR is polymorphic
grp = dplyr::group_by(tm3, gene)
tm4 = dplyr::summarise(grp, vartag = ifelse(sum(vartag=='Identical')==length(vartag), 'Identical','Polymorphic'))

# join the expression table with mapping table
ti3 = merge(ti2[,-c(2:5)], tm4, by='gene')
ti4 = gather(ti3[,c('gene','vartag','bpr.leaf','bpr.root')], tissue, bpr, -gene, -vartag)

p1 = ggplot(ti4) +
  #geom_histogram(aes(x = bpr, fill = vartag), position='dodge', alpha=0.5) + 
  geom_density(aes(x = bpr, fill = vartag), alpha=0.4, linetype='blank') + 
  geom_vline(xintercept = 0, color = 'forestgreen', size=0.1) +
  scale_x_continuous(name = 'log2(B/H)', limits=c(-5,5)) +
  #scale_y_discrete(name = 'RPKM') +
  scale_fill_brewer(name = 'Direction', palette = "Set1") + 
  facet_wrap(~tissue) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  #theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype = 0)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fo = sprintf("%s/11.bpr.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 5)