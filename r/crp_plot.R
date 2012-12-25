library(ggplot2)
library(data.table)
library(RColorBrewer)

ff = file.path(DIR_Misc2, "crp.annotation/family_info.tbl")
tf = read.table(ff, header=TRUE, sep="\t", as.is=T)
fams = unique(tf$cat)

org = "Mtruncatula"
org = "Mtruncatula_4.0"
org = "Athaliana"
fg = file.path(DIR_Misc3, "spada.crp", org, "31_model_SPADA", "61_final.tbl")
diro = file.path(DIR_Misc2, "crp.evo")

t1 = read.table(fg, header=TRUE, sep="\t", quote="", as.is=T)
colnames(t1)[which(colnames(t1)=="start")] = "beg"
t2 = t1[tolower(substr(t1$chr,1,3))=="chr" & tolower(t1$chr)!='chru', -which(colnames(t1)=="sequence")]

tli = t2
locs = as.numeric(substr(tli$chr, 4, 5)) * 1000000000 + (tli$beg + tli$end) / 2
names(locs) = tli$id
cl = locCluster(locs, 50000)
tlo = merge(tli, cl, by="id")

td = merge(tlo, tf, by.x="family", by.y="id")

pal <- colorRampPalette(brewer.pal(9,"Set1"))
p <- ggplot(data=td) +
	geom_point(mapping=aes(x=(beg+end)/2/1000000, y=cluster_y, colour=factor(cat, levels=fams)), shape=18, size=0.8) +
	scale_colour_manual(values=pal(length(unique(td$cat)))) +
	scale_x_continuous(name='Chr Position (Mbp)', expand=c(0.01, 0)) +
	scale_y_continuous(name='', expand=c(0.04, 0)) +
	facet_grid(chr ~ .) + 
	theme(legend.position='right', legend.title=element_blank()) +
	theme(axis.text.x=element_text(size=8, angle=0)) +
	theme(axis.text.y=element_blank(), axis.ticks=element_blank())
ggsave(p, filename=sprintf("%s/%s.png", diro, org), width=8, heigh=5)

write.table(g21[order(g21$pos),], file.path(dir, "02_cluster.tbl"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
