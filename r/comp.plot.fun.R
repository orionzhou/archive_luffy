## data processing functions
assign_block <- function(dfi, gap_prop=0.4, gap_len=5000) {
	if(ncol(dfi) == 3) { dfi = cbind(dfi, strand="+") }
	colnames(dfi)[1:4] = c('id', 'beg', 'end', 'strand')
	order.idx = order(dfi[,1], dfi[,2], dfi[,3])
	df = dfi[order.idx,]
	blocks = rep(0, nrow(df))
	cnt = 1
	for (i in 1:nrow(df)) {
		if( blocks[i] > 0 ) { next }
		blocks[i] = cnt
		idP=df[i,1]; begP=df[i,2]; endP=df[i,3]
		for (j in (i+1):nrow(df)) {
			id=df[j,1]; beg=df[j,2]; end=df[j,3]
			if( j == i | j > nrow(df) | id != idP ) { break }
			if( end < begP ) {
				gap = begP - end - 1
			} else if( endP < beg ) {
				gap = beg - endP + 1
			} else {
				gap = 0
			}
			blockBeg = min(beg, begP)
			blockEnd = max(end, endP)
			blockLen = blockEnd - blockBeg + 1
			if( gap/blockLen <= gap_prop | gap <= gap_len ) {
				blocks[j] = cnt
				idP=id; begP=blockBeg; endP=blockEnd
			} else {
				break
			}
		}
		cnt = cnt + 1
	}
	dfo = cbind(dfi, block_old=blocks[order(order.idx)])
	block_unique = unique(dfo$block_old)
	block_mapping = data.frame(block_old=block_unique, block=1:length(block_unique))
	dfo = cbind(dfo, idx=1:nrow(dfo))
	dfo = merge(dfo, block_mapping, by='block_old')
	dfo = dfo[order(dfo$idx),]
	dfo = dfo[, !colnames(dfo) %in% c('idx','block_old')]
	dfb = ddply(dfo, .(block), summarise, id=unique(id), beg=min(beg), end=max(end), strand=strand[which(end-beg == max(end-beg))[1]])
	list('df'=dfo, 'dfb'=dfb)
}

assign_block_mapping <- function(dfi, gap_len_q=10000, gap_len_h=100000) {
	order.idx = order(dfi[,1], dfi[,2], dfi[,3])
	df = dfi[order.idx,]
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
			} else if( abs(qBeg-qEndP)<gap_len_q & strd==strdP & (hId==hIdP) & ( (strd=="+" & abs(hBeg-hEndP)<gap_len_h) | (strd=="-" & abs(hEnd-hBegP)<gap_len_h) ) ) {
				blocks[j] = cnt
				qIdP=qId; qBegP=qBeg; qEndP=qEnd; strdP=strd; hIdP=hId; hBegP=hBeg; hEndP=hEnd;
			}
		}
		cnt = cnt + 1
	}
	dfo = cbind(dfi, block=blocks[order(order.idx)])
	dfb = ddply(dfo, .(block), summarise, 
		qId = unique(qId), qBeg = min(qBeg), qEnd = max(qEnd), strand = unique(strand),
		hId = unique(hId), hBeg = min(hBeg), hEnd = max(hEnd), qLen = max(qEnd)-min(qBeg)+1,
		hLen = max(hEnd)-min(hBeg)+1, qLen_aln = sum(qLen), hLen_aln = sum(hLen),
		pct = mean(pct), e = min(e), score = sum(score))
	list('df'=dfo, 'dfb'=dfb)
}

prepare_coord <- function(dfi) {
	df = dfi[order(dfi$block),]
	begr = as.integer( df$beg - 0.02 * (df$end - df$beg) )
	endr = as.integer( df$end + 0.02 * (df$end - df$beg) )
	df2 = cbind(df, begr=begr, endr=endr)
	df2$begr = apply(df2, 1, function(x) max(1,as.numeric(x['begr'])))
	df2$endr = apply(df2, 1, function(x) min(as.numeric(x['length']),as.numeric(x['endr'])))
	df3 = cbind(df2, len=df2$endr-df2$begr+1)

	len_total = sum(df3$len)
	itv = as.integer(0.01 * len_total)
	df5 = cbind(df3, beg.a = itv+1)
	if(nrow(df5) > 1) {
		for (i in 2:nrow(df5)) {
			df5$beg.a[i] = df5$beg.a[i-1] + df5$len[i-1] + itv + 1
		}
	}
	df6 = df5[,c('block','id','begr','endr','strand','len','beg.a')]
	colnames(df6)[3:4] = c('beg','end')
	df7 = cbind(df6, end.a=df6$beg.a+df6$len-1)
	df7
}

filter_assign_block <- function(dfi, dfb) {
	if(ncol(dfi) == 3) { dfi = cbind(dfi, strand="+") }
	colnames(dfi) = c('id','beg','end','strand')
	idxs.raw = c()
	blk.idxs = c()
	for (i in 1:nrow(dfb)) {
		id=dfb$id[i]; beg=dfb$beg[i]; end=dfb$end[i]; blk=dfb$block[i]
		idxss = which( dfi$id==id & ( (beg<=dfi$beg & dfi$beg<=end) | (beg<=dfi$end & dfi$end<=end) ) )
		if(length(idxss) > 0) {
			idxs.raw = c(idxs.raw, idxss)
			blk.idxs = c(blk.idxs, rep(i, length(idxss)))
		}
	}
	
	if(is.null(idxs.raw)) {idxs=NULL} else {idxs=sort(unique(idxs.raw))}
	blks = c()
	begs.a = c()
	ends.a = c()
	strands.a = c()
	for (i in idxs) {
		id=dfi[i,1]; beg=dfi[i,2]; end=dfi[i,3]; strand=dfi[i,4]
		
		blk.idxss = blk.idxs[ which(idxs.raw == i) ]
		blk.idx = blk.idxss[1]
		if( length(blk.idxss) > 1 ) {
			lens_ovlp = c()
			for (i in 1:length(blk.idxss) ) {
				len_ovlp = min(dfb$end[blk.idxss[i]], end) - max(dfb$beg[blk.idxss[i]], beg) + 1
				lens_ovlp = c(lens_ovlp, len_ovlp)
			}
			blk.idx = blk.idxss[ which(lens_ovlp == max(lens_ovlp))[1] ]
		}
		
		blk = dfb$block[blk.idx]
		blk.beg = dfb$beg[blk.idx]
		blk.end = dfb$end[blk.idx]
		blk.strand = dfb$strand[blk.idx]
		blk.beg.a = dfb$beg.a[blk.idx]
		blk.end.a = dfb$end.a[blk.idx]
		
		beg = max(beg, blk.beg)
		end = min(blk.end, end)
		beg.a = ifelse( blk.strand == "-", blk.end.a - (end - blk.beg), blk.beg.a + (beg - blk.beg) );
		end.a = ifelse( blk.strand == "-", blk.end.a - (beg - blk.beg), blk.beg.a + (end - blk.beg) );
		strand.a = ifelse(blk.strand == strand, "+", "-")
		
		blks = c(blks, blk)
		begs.a = c(begs.a, beg.a)
		ends.a = c(ends.a, end.a)
		strands.a = c(strands.a, strand.a)
	}
	data.frame(idx=idxs, block=blks, beg.a=begs.a, end.a=ends.a, strand.a=strands.a, stringsAsFactors=F)
}

get_ticks <- function(dfb, tick_itv) {
	dft = data.frame(block=c(),pos=c())	
	for (i in 1:nrow(dfb)) {
		tick_beg = pretty(c(dfb$beg[i], dfb$end[i]))[1]
		ticks = seq(tick_beg, dfb$end[i], by=tick_itv)
		if(length(ticks)==0) { ticks = tick_beg }
		dft = rbind(dft, data.frame(block=rep(dfb$block[i], length(ticks)), pos=ticks))
	}
	
	dft = merge(dft, dfb, by='block')
	dft = cbind(dft, pos.a=0)
	for (i in 1:nrow(dft)) {
		if(dft$strand[i] == "-") {
			dft$pos.a[i] = dft$end.a[i] - (dft$pos[i] - dft$beg[i])
		} else {
			dft$pos.a[i] = dft$beg.a[i] + (dft$pos[i] - dft$beg[i])
		}
	}
	dft = dft[dft$pos.a > 0, c('block','pos','strand','pos.a')]
	dfl = ddply(dft, .(block), summarise, beg.a=min(pos.a), end.a=max(pos.a))
	list(tick=dft, line=dfl)
}

data_preprocess <- function(tas, t_len, t_len_ref, t_gap, t_gap_ref, t_gen_ref) {
	dfm = assign_block_mapping(tas, 10000, 10000)$df
	dfm = dfm[order(dfm$hId, dfm$hBeg, dfm$hEnd), !colnames(dfm) %in% c('block')]

	list1 = assign_block(dfm[,c('qId', 'qBeg', 'qEnd', 'strand')])
	list2 = assign_block(dfm[,c('hId', 'hBeg', 'hEnd')])

	dfb1.1 = merge(list1$dfb, t_len, by='id')
	dfb1.2 = prepare_coord(dfb1.1)
	dfb1 = dfb1.2

	dfb2.1 = merge(list2$dfb, t_len_ref, by='id')
	dfb2.2 = prepare_coord(dfb2.1)
	dfb2 = dfb2.2

	max_len = max(dfb1$end.a+dfb1$beg.a[1]-1, dfb2$end.a+dfb2$beg.a[1]-1)
	tick_itv = diff( pretty(c(1,max(dfb1$len, dfb2$len)))[1:2] )
	xaxis1 = get_ticks(dfb1, tick_itv)
	xaxis2 = get_ticks(dfb2, tick_itv)

	dfm.a = filter_assign_block(dfm[,c('qId','qBeg','qEnd','strand')], dfb1)
	dfm.b = filter_assign_block(dfm[,c('hId','hBeg','hEnd')], dfb2)
	stopifnot(sort(dfm.a$idx) == 1:nrow(dfm), sort(dfm.b$idx)==1:nrow(dfm))
	dfm.2 = cbind(dfm, blk1=dfm.a$block, blk1.beg.a=dfm.a$beg.a, blk1.end.a=dfm.a$end.a, blk1.strand.a=dfm.a$strand.a, blk2=dfm.b$block, blk2.beg.a=dfm.b$beg.a, blk2.end.a=dfm.b$end.a, blk2.strand.a=dfm.b$strand.a, stringsAsFactors=F)
	comp = dfm.2

	annt1 = data.frame()

	df2 = t_gen_ref[,c('id','chr','beg','end','strand','note')]
	df2.1 = filter_assign_block(df2[,2:5], dfb2)
	df2.2 = cbind(df2[df2.1$idx,], df2.1[,-1])
	annt2 = df2.2

	dfg1.1 = filter_assign_block(t_gap[,1:3], dfb1)
	dfg1 = cbind(t_gap[dfg1.1$idx,], dfg1.1[,-1])

	dfg2.1 = filter_assign_block(t_gap_ref[,1:3], dfb2)
	dfg2 = cbind(t_gap_ref[dfg2.1$idx,], dfg2.1[,-1])

	color = c()
	fill = c('assembly gap'='black', 'gene(+)'='forestgreen', 'gene(-)'='dodgerblue')
	
	list(seg1=dfb1, seg2=dfb2, comp=comp, annt1=annt1, annt2=annt2, 
		xaxis1=xaxis1, xaxis2=xaxis2, gap1=dfg1, gap2=dfg2, 
		max_len=max_len, color=color, fill=fill)
}

## plotting functions
plot_segment <- function(df, y=unit(0.5,'npc'), col.p='red', col.n='blue', text.above=F, text.rot=0, vp=NULL) {
	colnames(df) = c('id','beg','end','strand')
	idxs.n = which(df$strand == "-")
	
	line.cols = rep(col.p, nrow(df))
	line.cols[idxs.n] = col.n
	grid.segments(
		x0 = unit(df$beg, 'native'), x1 = unit(df$end, 'native'),
		y0 = y, y1 = y, 
		gp = gpar(col=line.cols), vp = vp)
	
	for (i in 1:nrow(df)) {
		beg = df[i,2]; end = df[i,3]; strand = df[i,4]
		if(strand == '-') {
			x1 = unit(beg, 'native') + unit(5,'points')
			x2 = unit(beg, 'native')
		} else {
			x1 = unit(end, 'native') - unit(5,'points')
			x2 = unit(end, 'native')
		}
		if(i == 1) {
			arrows.x1 = x1
			arrows.x2 = x2
		} else {
			arrows.x1 = unit.c(arrows.x1, x1)
			arrows.x2 = unit.c(arrows.x2, x2)
		}
	}
	grid.segments(
		x0 = arrows.x1, x1 = arrows.x2,
		y0 = y - unit(3, 'points'), y1 = y, 
		gp = gpar(col=line.cols), vp = vp)
	grid.segments(
		x0 = arrows.x1, x1 = arrows.x2,
		y0 = y + unit(3, 'points'), y1 = y, 
		gp = gpar(col=line.cols), vp = vp)
	
	if(text.above) {
		text.y = y + unit(5, 'points')
	} else {
		text.y = y - unit(5, 'points')
	}
	text.offset = ifelse(text.above, unit(-10, "points"), unit(10, 'points'))
	grid.text( label = df$id, x = unit(df$beg, 'native'), 
		y = text.y, just = c("left","center"), 
		rot = text.rot, gp = gpar(cex=0.8, fontfamily="Helvetica"), vp = vp)
}
plot_xaxis <- function(xaxis, y=unit(0,'npc'), tick.above=F, vp=NULL) {
	dfl = xaxis$line
	dft = xaxis$tick
	grid.segments(
		x0 = unit(dfl$beg.a, 'native'),
		x1 = unit(dfl$end.a, 'native'),
		y0 = y, y1 = y, vp = vp)
	
	if( tick.above ) {
		tick.y = y + unit(3, 'points')
		text.y = y + unit(5, 'points')
		text.just = c('center', 'bottom')
	} else {
		tick.y = y - unit(3, 'points')
		text.y = y - unit(5, 'points')
		text.just = c('center', 'top')
	}
	grid.segments(
		x0 = unit(dft$pos.a, 'native'),
		x1 = unit(dft$pos.a, 'native'),
		y0 = y, y1 = tick.y, 
		vp = vp)
	grid.text( 
		label = dft$pos/1000,
		x = unit(dft$pos.a, 'native'), 
		y = text.y,	just = text.just, 
		gp = gpar(cex=0.7), vp = vp)
}
plot_feature_ds <- function(df, y=unit(0.5,'npc'), height=unit(5,'points'), fill.p='red', fill.n='blue', text.above=F, text.rot=0, vp=NULL) {
	colnames(df) = c('id','beg','end','strand')
	rect.fills = c()
	for (i in 1:nrow(df)) {
		if(df$strand[i] == "+") {
			rect.y = y + unit(2, 'points')
			rect.fill = fill.p
		} else {
			rect.y = y - height - unit(1, 'points')
			rect.fill = fill.n
		}
		if( i == 1 ) { rect.ys = rect.y } else { rect.ys = unit.c(rect.ys, rect.y) }
		rect.fills = c(rect.fills, rect.fill)
	}
	grid.rect( 
		x = unit(df$beg, 'native'), y = rect.ys,
		width = unit(df$end-df$beg, 'native'), 
		height = rep(height, nrow(df)), just = c('left', 'bottom'),
		gp = gpar(lwd=0, fill=rect.fills, alpha=0.9), vp = vp)	

	if(text.above) {
		text.y = y + unit(15, 'points')
	} else {
		text.y = y - unit(15, 'points')
	}
	grid.text( df$id, 
		x = unit(df$beg, 'native'), 
		y = text.y, just = c("left","center"), 
		rot = text.rot, gp = gpar(cex=0.8, fontfamily="Helvetica-Narrow"), vp = vp)
}
plot_feature <- function(df, y=unit(0.5,'npc'), height=unit(5,'points'), fill='grey', vp=NULL) {
	colnames(df) = c('id','beg','end')
	grid.rect( 
		x = unit(df$beg, 'native'), 
		y = y,
		width = unit(df$end-df$beg, 'native'), height=height,
		just=c('left','center'),
		gp = gpar(lwd=0, fill=fill, alpha=0.9), vp = vp)
}
plot_comparison <- function(comp, y1=unit(0.95, 'npc'), y2=unit(0.05, 'npc'), fill.p='skyblue1', fill.n='tomato', alpha=0.1, vp=NULL) {
	comp.xs = c()
	comp.fills = c()
	comp.ids = c()
	for (i in 1:nrow(comp) ) {
		if( comp$blk1.strand.a[i] == comp$blk2.strand.a[i] ) {
			comp.x = c(comp$blk1.beg.a[i], comp$blk1.end.a[i], comp$blk2.end.a[i], comp$blk2.beg.a[i])
			comp.fill = fill.p
			comp.id = rep(i, 4)
		} else {
			comp.x = c(comp$blk1.beg.a[i], comp$blk1.end.a[i], comp$blk2.beg.a[i], comp$blk2.end.a[i])
			comp.fill = fill.n
			comp.id = rep(i, 4)
		}
		comp.xs = c(comp.xs, comp.x)
		comp.fills = c(comp.fills, comp.fill)
		comp.ids = c(comp.ids, comp.id)
	}
	
	grid.polygon(
		x = unit(comp.xs, 'native'),
		y = rep(c(y1,y1,y2,y2), nrow(comp)),
		id = comp.ids,
		gp = gpar(fill=comp.fills, alpha=alpha, lwd=0),
		vp = vp)
}
plot_legend <- function(fill, x=unit(0.5,'npc'), y=unit(0.5,'npc'), height=unit(5,'points'), width=unit(30,'points'), vp=NULL) {
	n = length(fill)
	fill.labels = names(fill)
	grid.rect( 
		x = rep(x, n), 
		y = y - unit(seq(0, by=10, length.out=n), 'points'),
		width = width, height=height,
		just=c('left','center'),
		gp = gpar(lwd=0, fill=fill, alpha=0.9), vp = vp)
	grid.text( fill.labels, 
		x = x + width + unit(10, 'points'), 
		y = y - unit(seq(0, by=10, length.out=n), 'points'), 
		just = c("left","center"), 
		gp = gpar(cex=0.9, fontfamily="mono"), vp = vp)
}
plot_scale <- function(max_len, x=unit(0.5,'npc'), y=unit(0.5,'npc'), vp=NULL) {
	len = diff( pretty(1:max_len, 20)[1:2] )
	name = sprintf("%.0fkb", len/1000)
	grid.segments( 
		x0 = x, x1 = x + unit(len, 'native'),
		y0 = y, y1 = y,
		vp = vp)
	grid.segments( 
		x0 = unit.c(x, x + unit(len, 'native')),
		x1 = unit.c(x, x + unit(len, 'native')),
		y0 = rep(y, 2), 
		y1 = rep(y + unit(3, 'points'), 2), 
		vp = vp)
	grid.text( name,
		x = x + unit(len/2, 'native'), 
		y = y + unit(5, 'points'), 
		just = c("center","bottom"), 
		gp = gpar(cex=0.9, fontfamily="Helvetica"), vp = vp)
}

plot_final <- function(fn, dat, width=2000, height=1000, title="") {
	seg1=dat$seg1; seg2=dat$seg2; comp=dat$comp
	annt1=dat$annt1; annt2=dat$annt2
	xaxis1=dat$xaxis1; xaxis2=dat$xaxis2
	gap1=dat$gap1; gap2=dat$gap2
	color=dat$color; fill=dat$fill
	max_len=dat$max_len
	
	png(filename = fn, width = width, height = height)
	grid.newpage()
	vpf = viewport(	x = unit(0.5,"npc"), y = unit(0.5,"npc"), 
		width = unit(1, "npc") - unit(4, "lines"), height = unit(1, "npc") - unit(6,"lines"), 
		just=c("center", "center"))
	pushViewport(vpf)

	vp.layout = viewport(
		layout=grid.layout(3, 1, widths=unit(1, "npc"), heights=unit(c(0.3,0.4,0.3), "npc")))
	pushViewport(vp.layout)

	vp1 <- viewport(layout.pos.col=1, layout.pos.row=1, xscale=c(1,max_len), name='top')
	pushViewport(vp1); upViewport()
	vp2 <- viewport(layout.pos.col=1, layout.pos.row=3, xscale=c(1,max_len), name='bot')
	pushViewport(vp2); upViewport()
	vpm <- viewport(layout.pos.col=1, layout.pos.row=2, xscale=c(1,max_len), name='mapping')
	pushViewport(vpm); upViewport()

	# top 
	grid.rect(gp=gpar(fill='grey', alpha=0.2, lwd=0), vp = vp1)
	plot_segment(seg1[,c('id','beg.a','end.a','strand')], y=unit(10,'points'), text.above=T, text.rot=30, vp=vp1)
	plot_xaxis(xaxis1, y=unit(0,'npc'), tick.above=F, vp=vp1)

	# bottom
	grid.rect(gp=gpar(fill='grey', alpha=0.2, lwd=0), vp = vp2)
	plot_segment(seg2[,c('id','beg.a','end.a','strand')], y=unit(1,'npc')-unit(10,'points'), text.above=F, text.rot=-30, vp=vp2)
	plot_xaxis(xaxis2, y=unit(1,'npc'), tick.above=T, vp=vp2)

	if(nrow(annt2) > 0) {
		plot_feature_ds(annt2[,c('note','beg.a','end.a','strand.a')], y=unit(1,'npc')-unit(10,'points'), height=unit(5,'points'), fill.p=fill['gene(+)'], fill.n=fill['gene(-)'], text.above=F, text.rot=-40, vp=vp2)
	}

	# middle
	plot_comparison(comp, y1=unit(0.95, 'npc'), y2=unit(0.05, 'npc'), alpha=0.5, vp=vpm)
	if(nrow(gap1) > 0) {
		plot_feature(gap1[,c('id','beg.a','end.a')], y=unit(0.95,'npc')+unit(3, 'points'), height=unit(5,'points'), fill=fill['assembly gap'], vp=vpm)
	}
	if(nrow(gap2) > 0) {
		plot_feature(gap2[,c('id','beg.a','end.a')], y=unit(0.05,'npc')-unit(3, 'points'), height=unit(5,'points'), fill=fill['assembly gap'], vp=vpm)
	}

	# misc
	grid.text(title, 
		x = unit(0.5,'npc'), y = unit(1,'npc') - unit(2,'lines'), 
		gp = gpar(cex=1.5, fontface='bold', fontfamily='Helvetica'),
		vp=vpf)
	plot_legend(fill, x=unit(0.05,'npc'), y=unit(0.5,'npc'), vp=vp1)
	plot_scale(max_len, x=unit(0.9,'npc'), vp=vp1)

	dev.off()
}