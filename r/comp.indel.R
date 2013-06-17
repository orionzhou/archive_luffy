source("comp.plot.fun.R")

##### identify insertion / deletions
indel = data.frame()
gap1 = t_gap
gap2 = t_gap_ref
for (i in 2:nrow(t_aln)) {
	df = t_aln
	if( df$qId[i]==df$qId[i-1] & df$hId[i]==df$hId[i-1] & df$strand[i]==df$strand[i-1] ) {
		qId = df$qId[i]; hId = df$hId[i]; strand = df$strand[i]
		qBeg1 = df$qBeg[i-1]; qEnd1 = df$qEnd[i-1]; qBeg2 = df$qBeg[i]; qEnd2 = df$qEnd[i]
		hBeg1 = df$hBeg[i-1]; hEnd1 = df$hEnd[i-1]; hBeg2 = df$hBeg[i]; hEnd2 = df$hEnd[i]
		if( strand == "-" ) {
			qItv = qBeg2 - qEnd1 - 1
			hItv = hBeg1 - hEnd2 - 1
			qHasGap = sum( gap1$id==qId & gap1$beg >= qEnd1 & gap1$end <= qBeg2 )
			hHasGap = sum( gap2$id==hId & gap2$beg >= hEnd2 & gap2$end <= hBeg1 )
		} else {
			qItv = qBeg2 - qEnd1 - 1
			hItv = hBeg2 - hEnd1 - 1
			qHasGap = sum( gap1$id==qId & gap1$beg >= qEnd1 & gap1$end <= qBeg2 )
			hHasGap = sum( gap2$id==hId & gap2$beg >= hEnd1 & gap2$end <= hBeg2 )
		}
		type = ''
		if( abs(qItv) <= 10 & hItv >= 100 & hItv <= 10000 & hHasGap == 0) {
			type = 'del'
		} else if( abs(hItv) <= 10 & qItv >= 100 & qItv <= 10000 & qHasGap == 0) {
			type = "ins"
		}
		if(type != '') {
			indel.1 = data.frame('qId'=qId, 'qBeg1'=qBeg1, 'qEnd1'=qEnd1, 'qBeg2'=qBeg2, 'qEnd2'=qEnd2, 'strand'=strand, 
			'hId'=hId, 'hBeg1'=hBeg1, 'hEnd1'=hEnd1, 'hBeg2'=hBeg2, 'hEnd2'=hEnd2, 
			'qItv'=qItv, 'hItv'=hItv, 'qLen1'=df$qLen[i-1], 'qLen2'=df$qLen[i], 
			'hLen1'=df$hLen[i-1], 'hLen2'=df$hLen[i], 'pct1'=df$pct[i-1], 'pct2'=df$pct[i], 
			'score1'=df$score[i-1], 'score2'=df$score[i], 'type'=type)
			indel = rbind(indel, indel.1)
		}
	}
}
nrow(indel)
head(indel)

f31_01 = file.path(dir, "31_indel/01_raw.tbl")
write.table(indel, file=f31_01, col.names=T, row.names=F, sep="\t", quote=F)

##### read indels
dat1.name = "hm056"
dat1.dir = file.path("/home/youngn/zhoup/Data/misc3", dat1.name)

f31_01 = file.path(dat1.dir, "31_indel/01_raw.tbl")
indel = read.table(f31_01, header=TRUE, sep="\t", as.is=T)

idxs = c()
dfq = indel
tg = rbind(dat2$crp, dat2$nbs)
for( i in 1:nrow(dfq) ) {
	chr = dfq$hId[i];	beg = dfq$hEnd1[i]; end = dfq$hBeg2[i]
	flags = tg$chr==chr & ( (beg < tg$beg & tg$beg < end) | (beg < tg$end & tg$end < end) )
	if(sum(flags) > 0) {
		idxs = c(idxs, i)
	}
}
indel.1 = indel[idxs,]
nrow(indel.1)
head(indel.1)

dfp = indel.1
for (i in 1:nrow(dfp)) {
	chr = dfp$hId[i]
	beg = min( dfp$hBeg1[i], dfp$hBeg2[i] )
	end = max( dfp$hEnd1[i], dfp$hEnd2[i] )
	tas = aln[aln$hId==chr & ((beg<=aln$hBeg & aln$hBeg<=end) | (beg<=aln$hEnd & aln$hEnd<=end)), ]
	dat = data_preprocess(tas, dat1, dat2)
	fn = sprintf("%s/figs/indel_%d.png", dat1.dir, i)
	subtitle = sprintf("%s:%g-%gK", chr, beg/1000, end/1000)
	comp.plot(fn, dat, width=2000, height=1000, subtitle=subtitle)
}


idxs = c()
dir_pindel = "/home/youngn/zhoup/Data/misc3/hapmap_mt40/40_sv"
f_pindel_del = file.path(dir_pindel, "HM056_chr5_D.tbl")
td = read.table(f_pindel_del, head=T, sep="\t", as.is=T)
td.1 = td
td.2 = td.1[abs(td.1$svlen) >= 1 & abs(td.1$svlen) <= 10000, 1:7]

idxs = c()
dfq = td.2
tg = dat2$gene
tg = rbind(dat2$crp, dat2$nbs)
for( i in 1:nrow(dfq) ) {
	chr = dfq$chr[i]; beg = dfq$pos[i]; end =  dfq$end[i]
	flags = tg$chr==chr & ( (beg < tg$beg & tg$beg < end) | (beg < tg$end & tg$end < end) )
	if(sum(flags) > 0) { idxs = c(idxs, i) }
}
td.3 = dfq[idxs,]
nrow(td.3)
head(td.3)
