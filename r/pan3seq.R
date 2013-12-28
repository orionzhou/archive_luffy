library("GenomicRanges")
library("plyr")

# load query stats
org.q = "HM056"
dir.q = file.path('/home/youngn/zhoup/Data/genome', org.q)

fs.q = file.path(dir.q, '15_seqlen.tbl')
fg.q = file.path(dir.q, '16_gaploc.tbl')

ts.q = read.table(fs.q, sep='\t', header=T, as.is=T)
tg.q = read.table(fg.q, sep='\t', header=T, as.is=T)

si.q = Seqinfo(ts.q$id, seqlengths=ts.q$length)
gra.q = GRanges(seqnames=Rle(ts.q$id), ranges=IRanges(rep(1, nrow(ts.q)), end=ts.q$length), seqinfo=si.q)
grg.q = GRanges(seqnames=Rle(tg.q$id), ranges=IRanges(tg.q$beg, end=tg.q$end), seqinfo=si.q)

sum(width(gra.q))
sum(width(grg.q))

# load target stats
org.t = "HM101"
dir.t = file.path('/home/youngn/zhoup/Data/genome', org.t)

fs.t = file.path(dir.t, '15_seqlen.tbl')
fg.t = file.path(dir.t, '16_gaploc.tbl')

ts.t = read.table(fs.t, sep='\t', header=T, as.is=T)
tg.t = read.table(fg.t, sep='\t', header=T, as.is=T)

si.t = Seqinfo(ts.t$id, seqlengths=ts.t$length)
gra.t = GRanges(seqnames=Rle(ts.t$id), ranges=IRanges(rep(1, nrow(ts.t)), end=ts.t$length), seqinfo=si.t)
grg.t = GRanges(seqnames=Rle(tg.t$id), ranges=IRanges(tg.t$beg, end=tg.t$end), seqinfo=si.t)

sum(width(gra.t))
sum(width(grg.t))

# working dir
dir = sprintf('/home/youngn/zhoup/Data/misc3/%s_%s', org.q, org.t)

# generating novel sequences
pre="nov1"
fw = sprintf("%s/41_novelseq/%s.pre.gal", dir, pre)
tw = read.table(fw, sep='\t', header=T, as.is=T)[,1:17]
sum(tw$match) / (sum(tw$match) + sum(tw$misMatch))
fl = sprintf("%s/41_novelseq/%s.pre.gall", dir, pre)
tl = read.table(fl, sep='\t', header=T, as.is=T)

grm1 = GRanges(seqnames=Rle(tl$tId), ranges=IRanges(tl$tBeg, end=tl$tEnd), seqinfo=si.q)
grm2 = reduce(grm1)
sum(width(grm2))

grn1 = setdiff(gra.q, union(grm2, grg.q))
sum(width(grn1))

tn1 = data.frame(id=as.character(seqnames(grn1)), beg=as.numeric(start(grn1)), end=as.numeric(end(grn1)), len=as.numeric(width(grn1)))
tn2 = tn1[tn1$len>=1000,]
sum(tn2$len)

fn = sprintf("%s/41_novelseq/%s.tbl", dir, pre)
# write.table(tn2, fn, col.names=T, row.names=F, sep='\t', quote=F)
# seqextract.pl -i $data/genome/$org_q/11_genome.fa -o $pre.fa -n $pre.tbl


# plot novel segments length distribution
tmp1 = table(tn1$len)
tp.1 = data.frame(len=as.numeric(names(tmp1)), cnt=c(tmp1))
tp.2 = cbind(tp.1, sum=tp.1$len * tp.1$cnt)
tp.3 = tp.2[order(tp.2$len, decreasing=T),]
tp = cbind(tp.3, cumsum=cumsum(tp.3$sum))
plot(tp$len, tp$cumsum, type='l', xlim=c(0, 8000), main=org.q, xlab='segment length', ylab='cumsum of segments')
x=c(100,1000)
y=c(tp$cumsum[tp$len==100], tp$cumsum[tp$len==1000])
segments(0, y, x, y, col='blue')
segments(x, y, x, 0, col='blue')


# blastn validation
pre="nov1"
dirb=sprintf("%s/41_novelseq/%s.blast", dir, pre)
tvw = read.table(file.path(dirb, '12.gal'), sep='\t', header=T, as.is=T)
sum(tvw$match)/sum(tvw$match+tvw$misMatch)
tv = read.table(file.path(dirb, '12.gall'), sep='\t', header=T, as.is=T)
grv.1 = GRanges(seqnames=Rle(tv$qId), ranges=IRanges(tv$qBeg, end=tv$qEnd), seqinfo=si.q)
grv.2 = reduce(grv.1)
sum(width(grv.2))

# construct pseudo-chrs
# get scaffold order
tm = read.table(file.path(dir, "25_blat_final/35.gal"), header=T, sep="\t", as.is=T)[,1:17]
get_row_max_score <- function(df) { df[which.max(df[,'score']),] }
tm.s = ddply(tm, .(tId), get_row_max_score)
to.1 = data.frame(id=tm.s$tId, chr=tm.s$qId, pos=(tm.s$qBeg+tm.s$qEnd)/2, srd=tm.s$qSrd, idx=0, stringsAsFactors=F)
to = merge(to.1, ts.q, by='id', all=T)[,-6]
to = to[order(to$chr, to$pos, to$id),]
to$idx = 1:nrow(to)

tt = data.frame(id=to$id, cat_blat='mt', idx=to$idx, stringsAsFactors=F)
tt$cat_blat[is.na(to$chr)] = rep('unc', sum(is.na(to$chr)))

# blast NT
tb = read.table(file.path(dir, "51_pan3/15.tbl"), header=T, sep="\t")[,-c(18,19)]
tb.1 = ddply(tb, .(qId, cat), summarise, score=sum(score))
tb.2 = ddply(tb.1, .(qId), get_row_max_score)
tb = data.frame(id=tb.2$qId, cat_nt=as.character(tb.2$cat), stringsAsFactors=F)

# classify novelseq
ti = read.table(file.path(dir, "51_pan3/01.tbl"), header=T, sep="\t")[,-c(18,19)]
colnames(ti)[1] = 'chr'
makelocid <- function(x) sprintf("%s-%d-%d", x['chr'], as.numeric(x['beg']), as.numeric(x['end']))
ti = cbind(id=apply(ti, 1, makelocid), ti)

ti.1 = merge(ti, tb, by="id", all=T)
ti.2 = merge(ti.1, tt, by.x='chr', by.y='id', all.x=T)
ti.2$cat_nt[which(is.na(ti.2$cat_nt))] <- 'unc'
table(ti.2[,c('cat_nt', 'cat_blat')])
ti.3 = cbind(ti.2, cat='plant', stringsAsFactors=F)
ti.3$cat[which(ti.2$cat_blat=='unc' & ti.2$cat_nt=='foreign')] = 'unc'

tf1 = ti.3[ti.3$cat=='plant',]
tf1 = tf1[order(tf1$idx), c(1,3,4,5)]
cat(org.q, "specific:", nrow(tf1), "segments,", sum(tf1$len), "bp\n")
write.table(tf1, file.path(dir, "51_pan3/21_plant.tbl"), col.names=T, row.names=F, sep='\t', quote=F)
tf2 = ti.3[ti.3$cat=='unc',]
tf2 = tf2[order(tf2$idx), c(1,3,4,5)]
cat(org.q, "unc:", nrow(tf2), "segments,", sum(tf2$len), "bp\n")
write.table(tf2, file.path(dir, "51_pan3/22_unc.tbl"), col.names=T, row.names=F, sep='\t', quote=F)
