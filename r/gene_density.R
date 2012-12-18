a1 = readAssembly("mt_35")
a = droplevels(subset(a1, type!='centromere'))
c = droplevels(subset(a1, type=='centromere'))

f1 = file.path(DIR_Misc1, "window_stats/11_mt_35/51_stats.tbl")
d1 = read.table(f1, header=TRUE, sep="\t")
d1 = cbind(d1, position=(d1$beg+d1$end)/2, length=d1$end-d1$beg+1)
d2 = cbind(d1, percent_gene=d1$gene/d1$length, percent_TE=d1$transposable_element_gene/d1$length, dist=0)
for (chr in levels(d2$chr)) {
  idxCen = which(c$chr==chr)
  posCen = (c[idxCen, 'beg'] + c[idxCen, 'end']) / 2
  idxs = d2$chr==chr
  d2[idxs, 'dist'] = abs( (d2[idxs, 'beg'] + d2[idxs, 'end']) / 2 - posCen )
}

dirW = file.path(DIR_R, "mt_35")
write.table(d2, file.path(dirW, "01.tbl"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

p = ggplot(data = d2) + 
  geom_rect(mapping=aes(xmin=beg, xmax=end, ymin=0, ymax=percent_gene, fill='gene'), size=0) +
  geom_rect(mapping=aes(xmin=beg, xmax=end, ymin=percent_gene, ymax=percent_gene+percent_TE, fill='TE'), size=0) +
  layer(data=a, geom='rect', mapping=aes(xmin=beg, xmax=end, ymin=0.95, ymax=1, fill=type), geom_params=list(size=0)) +
  layer(data=c, geom='point', mapping=aes(x=(beg+end)/2, y=0.9, shape=type), geom_params=list(size=2)) +
  scale_fill_manual(values=c('phase 1'='aquamarine', 'phase 3'='aquamarine4', 'gene'='violet', 'TE'='tan'), legend=TRUE) +
  scale_shape_manual(values=c('centromere'=17), legend=TRUE) +
  labs(colour="Density", fill="Assembly", shape="") +
  scale_x_continuous(name='chr position (bp)', formatter='comma') +
  scale_y_continuous(name='', limits=c(0,1), formatter="percent") +
  opts(axis.text.y=theme_text(hjust=1, size=8, colour="blue")) +
  facet_grid(chr ~ .) 
ggsave(p, filename = file.path(dirW, "02_gd.png"), width=10, height=8)

p = ggplot(data = d1) + 
  geom_rect(mapping=aes(xmin=beg, xmax=end, ymin=0, ymax=cnt_TIR, fill='TIR'), size=0) +
  geom_rect(mapping=aes(xmin=beg, xmax=end, ymin=cnt_TIR, ymax=cnt_TIR+cnt_nonTIR, fill='nonTIR'), size=0) +
  geom_rect(mapping=aes(xmin=beg, xmax=end, ymin=cnt_TIR+cnt_nonTIR, ymax=cnt_TIR+cnt_nonTIR+cnt_CRP1040.1530, fill='NCR'), size=0) +
  geom_rect(mapping=aes(xmin=beg, xmax=end, ymin=cnt_TIR+cnt_nonTIR+cnt_CRP1040.1530, ymax=cnt_TIR+cnt_nonTIR+cnt_CRP1040.1530+cnt_CRP0000.1030, fill='nonNCR'), size=0) +
  scale_fill_brewer(palette='Set1') +	
  labs(fill="Gene Family") +
  scale_x_continuous(name='chr position (bp)', formatter='comma') +
  scale_y_continuous(name='Gene Number') +
  opts(axis.text.y=theme_text(hjust=1, size=8, colour="blue")) +
  facet_grid(chr ~ .) 
ggsave(p, filename = file.path(dirW, "03_gf.png"), width=10, height=8)

d1.1 = d1[,c(1:3,10:15)]
colnames(d1.1)[4:9] = c("gene", "TE", "non_nodule_DEFL", "nodule_DEFL", "nonTIR", "TIR")
d5 = melt(d1.1, id.vars=c('chr','beg','end'), variable_name='type')
for (i in 1:8) {
  chr = paste('chr', i, sep="")
  d5.1 = d5[d5$chr == chr,]
  p = ggplot(data = d5.1) +
    geom_rect(mapping=aes(xmin=beg, xmax=end, ymin=0, ymax=value, fill=type), size=0) + 
    scale_fill_brewer(palette='Dark2', legend=FALSE) + 
    scale_x_continuous(name='chr position (bp)', formatter='comma') +
    scale_y_continuous(name='Gene Number') +
    facet_grid(type ~ ., scales='free') +
    opts(axis.text.y=theme_text(hjust=1, size=7, colour="blue")) +
#  opts(strip.background = theme_rect(colour='NA',fill='NA')) +
    opts(strip.text.y = theme_text(size=7))
  ggsave(p, filename = sprintf("%s/03_gf_%s.png", dirW, chr), width=10, height=5)
}


