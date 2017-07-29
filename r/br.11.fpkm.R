require(dplyr)
require(GenomicRanges)

dirw = file.path(Sys.getenv("misc2"), "briggs")
dirw = '/home/springer/zhoux379/scratch/briggs2'

fi = file.path(dirw, '00.1.read.correct.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)[,1:5]

### get raw read count
gids = c()
for (i in 1:nrow(ti)) {
	sam = ti$SampleID[i]
	
	fh1 = sprintf("%s/32.htseq/%s.txt", dirw, sam)
	th1 = read.table(fh1, header = F, sep = "\t", as.is = T)
	ngene = nrow(th1) - 5
	gids0 = th1$V1[-c(ngene+1:5)]
	vals1 = th1$V2[-c(ngene+1:5)]
	
	fh2 = sprintf("%s/32.htseq/%s.as.txt", dirw, sam)
	th2 = read.table(fh2, header = F, sep = "\t", as.is = T)
	stopifnot(ngene == nrow(th2) - 5)
	stopifnot(gids0 == th2$V1[-c(ngene+1:5)])
	vals2 = th2$V2[-c(ngene+1:5)]
	if(i == 1) {
		gids = gids0
		to1 = data.frame(gid = gids, vals = vals1, stringsAsFactors = F)
		to2 = data.frame(gid = gids, vals = vals2, stringsAsFactors = F)
	} else {
		stopifnot(identical(gids, gids0))
		to1 = cbind(to1, vals = vals1)
		to2 = cbind(to2, vals = vals2)
	}
	colnames(to1)[i+1] = sam
	colnames(to2)[i+1] = sam
}

dim(to1)
fo1 = file.path(dirw, '32.rc.tsv')
write.table(to1, fo1, sep = "\t", row.names = F, col.names = T, quote = F)
fo2 = file.path(dirw, '32.rc.as.tsv')
write.table(to2, fo2, sep = "\t", row.names = F, col.names = T, quote = F)

### compute FPM
fi1 = file.path(dirw, '32.rc.tsv')
ti1 = read.table(fi1, header = T, sep = "\t", as.is = T)
total_reads1 = apply(ti1[,-1], 2, sum)
fi2 = file.path(dirw, '32.rc.as.tsv')
ti2 = read.table(fi2, header = T, sep = "\t", as.is = T)
total_reads2 = apply(ti2[,-1], 2, sum)
total_reads = total_reads1# + total_reads2

to1 = ti1
for (i in 1:length(total_reads)) {
	to1[,i+1] = to1[,i+1] / total_reads[i] * 1000000
}
fo1 = file.path(dirw, '33.fpm.tsv')
write.table(to1, fo1, sep = "\t", row.names = F, col.names = T, quote = F)

### compute FPKM
dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
f_gtb = file.path(dirg, "51.gtb")
f_tbl = file.path(dirg, "51.tbl")
tg = read.table(f_gtb, sep = "\t", header = T, as.is = T)[,1:2]
tt = read.table(f_tbl, sep = "\t", header = F, as.is = T)
colnames(tg) = c("tid", "gid")
colnames(tt) = c("chr", "beg", "end", "srd", "tid", "type", "fam")
tt2 = tt[tt$type %in% c('cds', 'utr5', 'utr3'),]
tg2 = merge(tg, tt2, by = 'tid')

gr = with(tg2, GRanges(seqnames = chr, ranges = IRanges(beg, end = end), gid = gid))
x = unlist(reduce(split(gr, elementMetadata(gr)$gid)))
tr = data.frame(gid = names(x), chr = seqnames(x), beg = start(x), end = end(x), stringsAsFactors = F)
grp = dplyr::group_by(tr, gid)
tr2 = dplyr::summarise(grp, len = sum(end - beg + 1))

to1 = merge(ti1, tr2, by = 'gid')
stopifnot(nrow(to1) == nrow(ti1))
for (i in 2:(ncol(to1)-1)) {
	to1[,i] = to1[,i] / (to1$len/1000)
}
dim(to1)
fo = file.path(dirw, '34.fpkm.tsv')
write.table(to1[,-ncol(to1)], fo, sep = "\t", row.names = F, col.names = T, quote = F)


### low count filter (obsolete)
fi1 = file.path(dirw, '34.fpkm.tsv')
fi2 = file.path(dirw, '34.fpkm.as.tsv')

ti1 = read.table(fi1, header = T, sep = "\t", as.is = T)
ti2 = read.table(fi2, header = T, sep = "\t", as.is = T)
e1 = asinh(ti1[,-1])
e2 = asinh(ti2[,-1])

nsam = ncol(e1)/3
lowexp = c()
for (i in 1:nsam) {
	n0 = apply(e1[,(3*(i-1)+1):(3*i)], 1, myfunc <- function(x) sum(x==0))
	e1s = e1[n0 == 2, (3*(i-1)+1):(3*i)]
	lowexp = c(lowexp, apply(e1s, 1, myfunc <- function(x) max(x)))
}

nsam = ncol(e1)/3
highexp = c()
for (i in 1:nsam) {
	n0 = apply(e1[,(3*(i-1)+1):(3*i)], 1, myfunc <- function(x) sum(x==0))
	e1s = e1[n0 == 0, (3*(i-1)+1):(3*i)]
	highexp = c(highexp, apply(e1s, 1, myfunc <- function(x) median(x)))
}

fo = file.path(dirw, 'stats/11.lowfilter.pdf')
pdf(fo, 6, 8)
par(mfrow=c(2,1))
hist(lowexp[lowexp<=6], 100, xlab = 'asinh(FPKM)', main = 'zero read count in 2 reps')
hist(highexp[highexp<=6], 100, xlab = 'asinh(FPKM)', main = 'non-zero read count in all 3 reps')
dev.off()

nsam = ncol(e2)/3
lowexp = c()
for (i in 1:nsam) {
	n0 = apply(e2[,(3*(i-1)+1):(3*i)], 1, myfunc <- function(x) sum(x==0))
	e2s = e2[n0 == 2, (3*(i-1)+1):(3*i)]
	lowexp = c(lowexp, apply(e2s, 1, myfunc <- function(x) max(x)))
}

highexp = c()
for (i in 1:nsam) {
	n0 = apply(e2[,(3*(i-1)+1):(3*i)], 1, myfunc <- function(x) sum(x==0))
	e2s = e2[n0 == 0, (3*(i-1)+1):(3*i)]
	highexp = c(highexp, apply(e2s, 1, myfunc <- function(x) median(x)))
}

fo = file.path(dirw, 'stats/11.lowfilter.as.pdf')
pdf(fo, 6, 8)
par(mfrow=c(2,1))
hist(lowexp[lowexp<=2], 100, xlab = 'asinh(FPKM)', main = 'zero read count in 2 reps')
hist(highexp[highexp<=2], 100, xlab = 'asinh(FPKM)', main = 'non-zero read count in all 3 reps')
dev.off()

### compute per-gene coverage (obsolete)
dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
f_gtb = file.path(dirg, "51.gtb")
f_tbl = file.path(dirg, "51.tbl")
tg = read.table(f_gtb, sep = "\t", header = T, as.is = T)[,1:2]
tt = read.table(f_tbl, sep = "\t", header = F, as.is = T)
colnames(tg) = c("tid", "gid")
colnames(tt) = c("chr", "beg", "end", "srd", "tid", "type", "fam")
tt2 = tt[tt$type %in% c('cds', 'utr5', 'utr3'),]
tg2 = merge(tg, tt2, by = 'tid')

gr = with(tg2, GRanges(seqnames = chr, ranges = IRanges(beg, end = end), gid = gid))
x = unlist(reduce(split(gr, elementMetadata(gr)$gid)))
tr = data.frame(gid = names(x), chr = seqnames(x), beg = start(x), end = end(x), stringsAsFactors = F)
grp = dplyr::group_by(tr, gid)
tr2 = dplyr::summarise(grp, len = sum(end - beg + 1))

fi = file.path(dirw, '32.rc.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

to2 = merge(ti, tr2, by = 'gid')
stopifnot(nrow(ti) == nrow(to2))
for (i in 2:(ncol(to2)-1)) {
	to2[,i] = to2[,i]*100 / to2$len
}
to3 = to2[,-ncol(to2)]
dim(to3)
fo = file.path(dirw, '32.cov.tsv')
write.table(to3, fo, sep = "\t", row.names = F, col.names = T, quote = F)

