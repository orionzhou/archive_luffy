library(ggplot2)
library(plyr)
library(RColorBrewer)

fg = file.path(DIR_Data, "genome/Mtruncatula_4.0/52_gap_loc.tbl")
tg = read.table(fg, header=TRUE, sep="\t", as.is=T)
#tg = tg[tg$len>=5000 & substr(tg$id,1,3)=="chr",]
colnames(tg)=c("hId", "hBeg", "hEnd", "len")

acc = "hm056"
#acc = "hm340"
dir = file.path(DIR_Misc3, acc)

f_len = file.path(dir, "11_seqlen.tbl")
t_len = read.table(f_len, header=TRUE, sep="\t", as.is=T)
colnames(t_len)=c("qId", "len_scaf")

f_gap = file.path(dir, "12_gaploc.tbl")
t_gap = read.table(f_gap, header=TRUE, sep="\t", as.is=T)
colnames(t_gap)=c("qId", "qBeg", "qEnd", "len")


dir = file.path(dir, '21_blastn')

f05 = file.path(dir, "05_tiled.tbl")
t05 = read.table(f05, header=TRUE, sep="\t", as.is=T)
t05  = t05[order(t05$qId, t05$qBeg),]


# identify insertion / deletions
find_indel <- function(i, df, tgQ, tgH) {
	if( df$qId[i]==df$qId[i-1] & df$hId[i]==df$hId[i-1] & df$strand[i]==df$strand[i-1] ) {
		qId = df$qId[i]
		hId = df$hId[i]
		strand = df$strand[i]
		if( t05$strand[i] == "-" ) {
			qItv = df$qBeg[i] - df$qEnd[i-1] - 1
			hItv = df$hBeg[i-1] - df$hEnd[i] - 1
			qHasGap = sum( tgQ$qId==qId & tgQ$qBeg >= df$qEnd[i-1] & tgQ$qEnd <= df$qBeg[i] )
			hHasGap = sum( tgH$hId==hId & tgH$hBeg >= df$hEnd[i-1] & tgH$hEnd <= df$HBeg[i] )
		} else {
			qItv = df$qBeg[i] - df$qEnd[i-1] - 1
			hItv = df$hBeg[i] - df$hEnd[i-1] - 1
			qHasGap = sum( tgQ$qId==qId & tgQ$qBeg >= df$qEnd[i-1] & tgQ$qEnd <= df$qBeg[i] )
			hHasGap = sum( tgH$hId==hId & tgH$hBeg >= df$hEnd[i] & tgH$hEnd <= df$HBeg[i-1] )
		}
		type = ''
		if( qItv < 10 & hItv > 100 & hHasGap == 0) {
			type = 'del'
		} else if( qItv > 100 & hItv < 10 & qHasGap == 0) {
			type = "ins"
		}
		if(type != '') {
			c('qId'=qId, 'qBeg1'=df$qBeg[i-1], 'qEnd1'=df$qEnd[i-1], 'qBeg2'=df$qBeg[i], 'qEnd2'=df$qEnd[i], 'strand'=strand, 'hId'=hId, 'hBeg1'=df$hBeg[i-1], 'hEnd1'=df$hEnd[i-1], 'hBeg2'=df$hBeg[i], 'hEnd2'=df$hEnd[i], 'qItv'=qItv, 'hItv'=hItv, 'qLen1'=df$qLen[i-1], 'qLen2'=df$qLen[i], 'hLen1'=df$hLen[i-1], 'hLen2'=df$hLen[i], 'pct1'=df$pct[i-1], 'pct2'=df$pct[i], 'score1'=df$score[i-1], 'score2'=df$score[i], 'type'=type)
		}
	}
}
t11 = ldply(2:nrow(t05), find_indel, t05, t_gap, tg)
nrow(t11)
head(t11)

dir = file.path(dir, "../31_indel")
f01 = file.path(dir, "01_raw.tbl")
write.table(t11, file=f01, col.names=T, row.names=F, sep="\t", quote=F)

t01 = read.table(f01, header=TRUE, sep="\t", as.is=T)

f_crp = file.path(DIR_Misc4, "spada.crp.Mtruncatula_4.0/31_model_evaluation/61_final.gtb")
t_crp = read.table(f_crp, header=TRUE, sep="\t", as.is=T)[,1:6]

overlap_indel_gene <- function(i, ti, tg) {
	chr = ti$hId[i]
	beg = ti$hEnd1[i]
	end = ti$hBeg2[i]
	flags = tg$chr==chr & ( (beg < tg$beg & tg$beg < end) | (beg < tg$end & tg$end < end) )
	if(sum(flags) > 0) {
		c('chr'=chr, 'beg'=beg, 'end'=end, 'size'=ti$hItv[i], 'gene'=tg$id[flags][1])
	}
}
to = ldply(1:nrow(t01), overlap_indel_gene, t01, t_crp)
to2 = to[as.numeric(to$hItv)>0 & as.numeric(to$hItv)<10000,]
head(to2)

# assign alignment blocks
assign_block <- function(df, qGap=10000, hGap=100000) {
	blocks = rep(0, nrow(df))
	cnt = 1
	for (i in 1:nrow(df)) {
		if( blocks[i] > 0 ) { next }
		blocks[i] = cnt
		qIdP=df$qId[i]; qBegP=df$qBeg[i]; qEndP=df$qEnd[i]; strdP=df$strand[i];
		hIdP=df$hId[i]; hBegP=df$hBeg[i]; hEndP=df$hEnd[i]
		for (j in (i+1):nrow(df)) {
			qId=df$qId[j]; qBeg=df$qBeg[j]; qEnd=df$qEnd[j]; strd=df$strand[j];
			hId=df$hId[j]; hBeg=df$hBeg[j]; hEnd=df$hEnd[j]
			if( j == i | j > nrow(df) ) { next }
			if( qId != qIdP ) {
				break
			} else if( abs(qBeg-qEndP)<qGap & strd==strdP & (hId==hIdP) & ( (strd=="+" & abs(hBeg-hEndP)<hGap) | (strd=="-" & abs(hEnd-hBegP)<hGap) ) ) {
				blocks[j] = cnt
				qIdP=qId; qBegP=qBeg; qEndP=qEnd; strdP=strd; hIdP=hId; hBegP=hBeg; hEndP=hEnd;
			}
		}
		cnt = cnt + 1
	}
	blocks
}
blocks = assign_block(t05)
t11 = cbind(t05, block=blocks)
f11 = file.path(dir, "11_blocks.tbl")
write.table(t11, file=f11, col.names=T, row.names=F, sep="\t", quote=F)

t12 = ddply(t11, .(qId, block), summarise, 
	qBeg = min(qBeg),
	qEnd = max(qEnd),
	strand = unique(strand),
	hId = unique(hId),
	hBeg = min(hBeg),
	hEnd = max(hEnd),
	qLen = max(qEnd)-min(qBeg)+1,
	hLen = max(hEnd)-min(hBeg)+1,
	qLen_aln = sum(qLen),
	hLen_aln = sum(hLen),
	pct = mean(pct),
	e = min(e),
	score = sum(score)
)
t12b = cbind(t12, qGap=t12$qLen-t12$qLen_aln, hGap=t12$hLen-t12$hLen_aln)
t12c = t12b[t12b$hGap > 100000,]
f12 = file.path(dir, "12_info_block.tbl")
write.table(t12b, file=f12, col.names=T, row.names=F, sep="\t", quote=F)

# get block information
f11 = file.path(dir, "11_blocks.tbl")
t11 = read.table(f11, header=TRUE, sep="\t", as.is=T)

tw2 = merge(tw, t_len, by='qId')
tw2 = cbind(tw2, pct_cov=tw2$qLen_aln/tw2$len_scaf)

sum_tw <- function(cutoff, df) {
	df = df[df$pct_cov >= cutoff,]
	c('num_qId'=length(unique(df$qId)), 'qLen_aln'=sum(df$qLen_aln), 'hLen_aln'=sum(df$hLen_aln), 'qLen'=sum(df$qLen))
}
ldply(seq(0,0.7,0.1), sum_tw, tw2)

tw3 = tw2[tw2$pct_cov >= 0.1,]
tws = ddply(tw3, .(qId), summarise, qLen_aln=sum(qLen_aln), hLen_aln=sum(hLen_aln), n_block=length(hId), len_scaf=unique(len_scaf), pct_cov=sum(pct_cov))


fw = file.path(dir, "09_wide.tbl")
tw = read.table(fw, header=TRUE, sep="\t", as.is=T)

fl = file.path(dir, "08_long.tbl")
tl = read.table(fl, header=TRUE, sep="\t", as.is=T)

fn = file.path(dir, "36_desc.tbl")
tn = read.table(fn, header=TRUE, sep="\t", as.is=T)



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

# plot genome distribution
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
has_gap_q = function(row, tgq) nrow( tgq[tgq$qId==row[1] & tgq$qBeg>as.numeric(row[2]) & tgq$qEnd<as.numeric(row[3]),] ) > 0
tws = tw$qId[ apply(tw, 1, has_gap, tg) ]

plot_scaffold <- function(id, tl, tw, tg, tqp, dir) {
	tws = tw[tw$qId==id,]
	tls = tl[tl$qId==id,]
	tls_f = tls[tls$strand=="+",]
	tls_r = tls[tls$strand=="-",]
	tgs = rbind.fill( apply(tws, 1, function(row, tg) tg[tg$hId==row[5] & tg$hBeg>as.numeric(row[6]) & tg$hEnd<as.numeric(row[7]), ], tg) )
	tgqs = rbind.fill( apply(tws, 1, function(row, tgq) tgq[tgq$qId==row[1] & tgq$qBeg>as.numeric(row[2]) & tgq$qEnd<as.numeric(row[3]), ], tgq) )
	p <- ggplot()
	if(nrow(tls_f) > 0) {
		p <- p + 
		geom_segment(data=tls[tls$strand=="+",], mapping=aes(x=hBeg/1000000, xend=hEnd/1000000, y=qBeg/1000000, yend=qEnd/1000000, color=pct), size=0.8)
	}
	if(nrow(tls_r) > 0) {
		p <- p +
		geom_segment(data=tls[tls$strand=="-",], mapping=aes(x=hBeg/1000000, xend=hEnd/1000000, y=qEnd/1000000, yend=qBeg/1000000, color=pct), size=0.8)
	}
	p <- p +
		scale_x_continuous(name='Position (Mbp) on Chr', expand=c(0.01, 0)) +
		scale_y_continuous(name='Position (Mbp) on Scaffold', expand=c(0.01, 0)) +
		scale_colour_continuous() + 
		facet_grid(qId ~ hId) + 
		theme_bw() +
		theme(legend.position='right') + 
		labs(colour="Percent Identity", fill="")
	if(length(tgs) > 0) {	
		p <- p +
		layer(data=tgs, geom='rect', mapping=aes(xmin=hBeg/1000000, xmax=hEnd/1000000, ymin=-Inf, ymax=Inf, fill="Mt4.0 Gap"), alpha=0.4)
	}
	if(length(tgqs) > 0) {	
		p <- p +
		layer(data=tgqs, geom='rect', mapping=aes(xmin=-Inf, xmax=Inf, ymin=qBeg/1000000, ymax=qEnd/1000000, fill="Scaffold Gap"), alpha=0.4)
	}
	p <- p +
		scale_fill_manual(values=c("red", "blue"))
	ggsave(p, filename=sprintf("%s/../figs/%s.png", dir, id), width=6.5, height=5)
}
nums = c(0,1,119,189,359,432,1000,1180,1181,1186,1188,1191,1203,1824,3189,4144)
for (num in nums) {
	id = sprintf("scaffold_%d", num)
	plot_scaffold(id, tl, tw, tg, tgq, dir)
}

chr = 'chr7'
beg = 20022555
end = 20116342
plot_region(chr, beg, end, tl, tw, tg, tgq, dir)
plot_region <- function(chr, beg, end, tl, tw, tg, tgq, dir) {
	tws = tw[tw$hId==chr & ( (beg<=tw$hBeg & tw$hBeg<=end) | (beg<=tw$hEnd & tw$hEnd<=end) ),]
	tls = tl[tl$hId==chr & ( (beg<=tl$hBeg & tl$hBeg<=end) | (beg<=tl$hEnd & tl$hEnd<=end) ),]
	tls_f = tls[tls$strand=="+",]
	tls_r = tls[tls$strand=="-",]
	tgs = rbind.fill( apply(tws, 1, function(row, tg) tg[tg$hId==row[5] & tg$hBeg>as.numeric(row[6]) & tg$hEnd<as.numeric(row[7]), ], tg) )
	tgqs = rbind.fill( apply(tws, 1, function(row, tgq) tgq[tgq$qId==row[1] & tgq$qBeg>as.numeric(row[2]) & tgq$qEnd<as.numeric(row[3]), ], tgq) )
	p <- ggplot()
	if(nrow(tls_f) > 0) {
		p <- p + 
		geom_segment(data=tls[tls$strand=="+",], mapping=aes(x=hBeg/1000000, xend=hEnd/1000000, y=qBeg/1000000, yend=qEnd/1000000, color=pct), size=0.8)
	}
	if(nrow(tls_r) > 0) {
		p <- p +
		geom_segment(data=tls[tls$strand=="-",], mapping=aes(x=hBeg/1000000, xend=hEnd/1000000, y=qEnd/1000000, yend=qBeg/1000000, color=pct), size=0.8)
	}
	p <- p +
		scale_x_continuous(name='Position (Mbp) on Chr', expand=c(0.01, 0)) +
		scale_y_continuous(name='Position (Mbp) on Scaffold', expand=c(0.01, 0)) +
		scale_colour_continuous() + 
		facet_grid(qId ~ hId) + 
		theme_bw() +
		theme(legend.position='right') + 
		labs(colour="Percent Identity", fill="")
	if(length(tgs) > 0) {	
		p <- p +
		layer(data=tgs, geom='rect', mapping=aes(xmin=hBeg/1000000, xmax=hEnd/1000000, ymin=-Inf, ymax=Inf, fill="Mt4.0 Gap"), alpha=0.4)
	}
	if(length(tgqs) > 0) {	
		p <- p +
		layer(data=tgqs, geom='rect', mapping=aes(xmin=-Inf, xmax=Inf, ymin=qBeg/1000000, ymax=qEnd/1000000, fill="Scaffold Gap"), alpha=0.4)
	}
	p <- p +
		scale_fill_manual(values=c("red", "blue"))
	ggsave(p, filename=sprintf("%s/../figs/%s_%d_%d.png", dir, chr, beg, end), width=8, height=5)
}