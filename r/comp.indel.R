##### identify insertion / deletions
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
t11 = ldply(2:nrow(t05), find_indel, t05, t_gap, t_gap_ref)
nrow(t11)
head(t11)

f01 = file.path(dir, "31_indel/01_raw.tbl")
write.table(t11, file=f01, col.names=T, row.names=F, sep="\t", quote=F)

t01 = read.table(f01, header=TRUE, sep="\t", as.is=T)

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
to2 = to[as.numeric(to$size)>0 & as.numeric(to$size)<10000,]
head(to2)