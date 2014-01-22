require(ggplot2)
require(Hmisc)
dir = file.path(DIR_Misc2, "genefam")

ta1 = read.table(file.path(dir, "31_merged.gtb"), header=T, as.is=T, sep="\t", quote="")

ta2 = cbind(ta1[,c('id','chr','beg','end','strand')], cat=ta1$cat1)
ta2$cat = as.character(ta2$cat)
idxs_cat2 = ta1$cat2!=''
ta2$cat[idxs_cat2] = ta1$cat2[idxs_cat2]

table(ta2$cat)
write.table(ta2, file.path(dir, "41_genefam.tbl"), sep="\t", col.names=T, row.names=F, quote=F)

ta1 = read.table(file.path(dir, "41_genefam.tbl"), sep="\t", header=T, as.is=T, quote="")
a1 = read.table(file.path(dir, "51_alignability.tbl"), sep="\t", header=T, as.is=T)
ta2 = merge(ta1, a1, by.x="id", by.y="id")
ta3 = cbind(ta2, pctC=ta2$lenC_a/ta2$lenC, pctM=ta2$lenM_a/ta2$lenM)

ss = function(df) { 
  x1 = df$pctC
  x2 = df$pctM
  q1 = quantile(x1, c(0.25, 0.5, 0.75))
  names(q1) = c('cds.0.25', 'cds.0.5', 'cds.0.75')
  q2 = quantile(x2, c(0.25, 0.5, 0.75))
  names(q2) = c('gene.0.25', 'gene.0.5', 'gene.0.75')
  c(sample.size=length(x1), cds.mean=mean(x1), cds.sd=sd(x1), q1, gene.mean=mean(x2), gene.sd=sd(x2), q2)
}
ddply(ta3, .(cat), ss)

p = ggplot(ta3, aes(x=cat, y=pctC, fill=cat)) +
#  stat_sum_df("mean_sdl", mult=1) + 
  geom_boxplot() + 
  scale_fill_brewer(palette='Set1') +  
  labs(fill="Gene Family") +
  scale_y_continuous(name='% alignability') +
  scale_x_discrete(name='') +
  theme(axis.text.x=theme_text(size=8, colour="blue", angle=10))
ggsave(p, filename = file.path(dir, "52_alignability.png"), width=8, height=6)


