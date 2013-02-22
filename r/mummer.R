library(ggplot2)
library(plyr)
library(data.table)
library(RColorBrewer)

fg = file.path(DIR_Data, "genome/Mtruncatula_4.0/51_stat.tbl")
tg = read.table(fg, header=TRUE, sep="\t", as.is=T)
#tg = tg[tg$len>=5000 & substr(tg$id,1,3)=="chr",]
colnames(tg)=c("hId", "hBeg", "hEnd", "len", "type")

dir = file.path(DIR_Misc3, "hm056/mummer")

# read data
fs = file.path(dir, "../hm056_seqlen.tbl")
ts = read.table(fs, header=TRUE, sep="\t", as.is=T)

ff = file.path(dir, "06.tbl")
tf = read.table(ff, header=TRUE, sep="\t", as.is=T)

fw = file.path(dir, "09_wide.tbl")
tw = read.table(fw, header=TRUE, sep="\t", as.is=T)

fl = file.path(dir, "08_long.tbl")
tl = read.table(fl, header=TRUE, sep="\t", as.is=T)

ts = cbind(ts, status=0)
id_mapped = unique(tf$qId)
ts$status[which(ts$id %in% id_mapped)] = 1
id_mapped_good = unique(tw$qId)
ts$status[which(ts$id %in% id_mapped_good)] = 2
write.table(ts, file.path(dir, "11_status.tbl"), col.names=T, row.names=F, sep="\t", quote=F)

scaffold_stat <- function(ids, ts) {
	lens = ts$length[ts$id %in% ids]
	c('n'=length(ids), 'total'=sum(lens), 'mean'=mean(lens), 'median'=median(lens))
}
sapply(list(all=id_all, mummer=, unmapped=id_unmapped, mummerG=id_mummerG, mummerB=id_mummerB), scaffold_stat, ts)

# plot raw scaffold mapping
plot_scaffold_raw <- function(id, tf, ts, dir) {
	tfs = tf[tf$qId==id,]
	lenS = ts$length[ts$id==id]
	p <- ggplot(data=tfs) +
		geom_segment(mapping=aes(x=qBeg/1000000, xend=qEnd/1000000, y=hBeg/1000000, yend=hEnd/1000000, color=pct_idty), size=0.5) +
		scale_x_continuous(name='Position (Mbp) on Scaffold', expand=c(0.01, 0), limits=c(0,lenS/1000000)) +
		scale_y_continuous(name='Position (Mbp) on Chr', expand=c(0.01, 0)) +
		scale_colour_continuous() + 
		facet_grid(hId ~ qId) + 
		theme_bw() +
		theme(legend.position='right') + 
		labs(colour="Percent Identity", fill="")
	ggsave(p, filename=sprintf("%s/figs/raw_%s.png", dir, id), width=6, height=8)
}
for (num in c(1:10,101:110)) {
	id = sprintf("scaffold_%d", num)
	id = id_mummerB[num]
	plot_scaffold_raw(id, tf, ts, dir)
}

# plot scaffold size distribution
tmp = cut(ts$length/1000, breaks=c(0,1,5,10,50,100,500,1000,5000))
p = ggplot(data.frame(size=tmp)) +
	geom_bar(aes(x=factor(size)), width=0.7) + 
	scale_x_discrete(name="Scaffold Size (kb)") + 
	scale_y_continuous(name="") +
	theme(axis.text.x = element_text(angle=15, size=8))
ggsave(file.path(dir, "figs/01_scaffold_size.png"), p, width=5, height=4)

# plot genome distribution of mummer mapping
p <- ggplot(data=tw) +
	geom_rect(mapping=aes(xmin=hBeg/1000000, xmax=hEnd/1000000, ymin=0, ymax=1, fill=hId)) +  
	layer(data=tg, geom='rect', mapping=aes(xmin=hBeg/1000000, xmax=hEnd/1000000, ymin=-1, ymax=0), geom_params=list()) +
	scale_x_continuous(name='Chr Position (Mbp)', expand=c(0.01, 0)) +
	scale_y_continuous(name='', expand=c(0.04, 0)) +
	facet_grid(hId ~ .) + 
	theme(legend.position='right', legend.title=element_blank()) +
	theme(axis.text.x=element_text(size=8, angle=0)) +
	theme(axis.text.y=element_blank(), axis.ticks=element_blank())
ggsave(p, filename=file.path(dir, "figs/03_coverage.png"), width=7, height=5)

# plot scaffold mapping with Mt4.0 gap
has_gap = function(row, tg) nrow( tg[tg$hId==row[5] & tg$hBeg>as.numeric(row[6]) & tg$hEnd<as.numeric(row[7]),] ) > 0
tws = tw$qId[ apply(tw, 1, has_gap, tg) ]

plot_scaffold <- function(id, tl, tw, tg, ts, dir) {
	tws = tw[tw$qId==id,]
	tls = tl[tl$qId==id,]
	tls_f = tls[tls$strand=="+",]
	tls_r = tls[tls$strand=="-",]
	tgs = rbind.fill( apply(tws, 1, function(row, tg) tg[tg$hId==row[5] & tg$hBeg>as.numeric(row[6]) & tg$hEnd<as.numeric(row[7]), ], tg) )
	lenS = ts$length[ts$id==id]
	p <- ggplot()
	if(nrow(tls_f) > 0) {
		p <- p + 
		geom_segment(data=tls[tls$strand=="+",], mapping=aes(x=qBeg/1000000, xend=qEnd/1000000, y=hBeg/1000000, yend=hEnd/1000000, color=pct_idty), size=0.8)
	}
	if(nrow(tls_r) > 0) {
		p <- p +
		geom_segment(data=tls[tls$strand=="-",], mapping=aes(x=qBeg/1000000, xend=qEnd/1000000, y=hEnd/1000000, yend=hBeg/1000000, color=pct_idty), size=0.8)
	}
	p <- p +
		scale_x_continuous(name='Position (Mbp) on Scaffold', expand=c(0.01, 0), limits=c(0,lenS/1000000)) +
		scale_y_continuous(name='Position (Mbp) on Chr', expand=c(0.01, 0)) +
		scale_colour_continuous() + 
		facet_grid(hId ~ qId, scales="free") + 
		theme_bw() +
		theme(legend.position='right') + 
		labs(colour="Percent Identity", fill="")
	if(length(tgs) > 0) {	
		p <- p +
		layer(data=tgs, geom='rect', mapping=aes(xmin=0, ymin=hBeg/1000000, ymax=hEnd/1000000, fill="Mt4.0 Gap"), xmax=lenS/1000000, alpha=0.4) + 
		scale_fill_manual(values=c("red"))
	}
	ggsave(p, filename=sprintf("%s/figs/%s.png", dir, id), width=6.5, height=5)
}
for (num in c(0,1,119,189,359,432,1000,1180,1181,1186,1188,1191,1203,1824,3189,4144)) {
	id = sprintf("scaffold_%d", num)
	plot_scaffold(id, tl, tw, tg, ts, dir)
}