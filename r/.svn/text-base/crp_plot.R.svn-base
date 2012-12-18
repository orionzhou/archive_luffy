dir = file.path(DIR_Misc2, "crp_plot")
fi = file.path(DIR_Misc2, "genefam/41_genefam.tbl")
g1 = read.table(fi, header=TRUE, sep="\t", quote="", as.is=T)
g2 = g1[g1$type1=='crp_gene',]

g3 = g2[as.character(g2$chr)!="chrU",]

freq = table(g3$note)
fams = names(freq[freq>20])
g4 = cbind(g3, family=g3$note)
g4$family=as.character(g4$family)
idxs_other = ! g4$family %in% fams
g4$family[idxs_other] = rep("other", sum(idxs_other))

g10 = g4

pos = as.numeric(substr(g10$chr, 4, 5)) * 1000000000 + (g10$beg + g10$end) / 2
names(pos) = g10$id
cl = locCluster(pos, 100000)

g20 = merge(g10, cl, by="id")

p <- ggplot(data=g20) +
	geom_point(mapping=aes(x=(beg+end)/2, y=cluster_y, colour=family), shape=18, size=0.7) +
	scale_colour_brewer(palette='Set1') +
	scale_x_continuous(name='chr position (bp)', formatter='comma', limits=c(1,48000000)) +
	scale_y_continuous(name='', breaks=seq(0,22,10), limits=c(1,22)) +
	facet_grid(chr ~ .) + 
	opts(legend.title=theme_blank()) +
	opts(axis.text.x=theme_text(size=8, angle=0)) +
	opts(axis.text.y=theme_blank(), axis.ticks=theme_blank(), axis.ticks.margin=unit(0,'cm'))
ggsave(p, filename=file.path(dir, '01.png'), width=8, heigh=5)

g21 = merge(g2, cl, by="id", all.x=all)
write.table(g21[order(g21$pos),], file.path(dir, "02_cluster.tbl"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

