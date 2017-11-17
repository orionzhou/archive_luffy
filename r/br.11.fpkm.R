require(dplyr)
require(GenomicRanges)

dirw = '/home/springer/zhoux379/scratch/briggs'
dirw = file.path(Sys.getenv("misc2"), "briggs")

fh = file.path(dirw, '00.1.read.correct.tsv')
th = read.table(fi, header = T, sep = "\t", as.is = T)[,1:5]

### get raw read count
gids = c()
for (i in 1:nrow(th)) {
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

### compute sample-level FPM
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

### compute sample-level FPKM
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


### tissue-level FPM & FPKM
sids_rm = c("BR045", "BR032", "BR026", "BR039", "BR042", "BR083", "BR095")
th = th[! th$SampleID %in% sids_rm,]

fr = file.path(dirw, '32.rc.tsv')
tr = read.table(fr, header = T, sep = "\t", as.is = T)
trl = gather(tr, sid, rc, -gid)
trl2 = merge(trl, th[,c("SampleID", "Tissue", "Genotype")], by.x = 'sid', by.y = 'SampleID')
grp = dplyr::group_by(trl2, gid, Tissue, Genotype)
trl3 = dplyr::summarise(grp, rc = sum(rc))
head(trl3)
nrow(trl3)

grp = dplyr::group_by(trl3, Tissue, Genotype)
trl4 = dplyr::summarise(grp, trc = sum(rc))
head(trl4)
nrow(trl4)

trl5 = merge(trl3, trl4, by = c("Tissue", "Genotype"))
trl6 = within(data.frame(trl5), {
    fpm = rc/trc*1000000
})
head(trl6)

### FPKM
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

to = trl6[,c('gid','Tissue','Genotype','fpm')]
to2 = merge(to, tr2, by = 'gid')
stopifnot(nrow(to) == nrow(to2))
to3 = cbind(to2, fpkm = to2$fpm / (to2$len / 1000))
to4 = to3[,-5]
dim(to4)
fo = file.path(dirw, '35.long.tsv')
write.table(to4, fo, sep = "\t", row.names = F, col.names = T, quote = F)


### filtering
fi = file.path(dirw, "35.long.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

grp = dplyr::group_by(ti, gid, Genotype)
ti2 = dplyr::summarise(grp, nexp = sum(fpkm>=1))
ti3 = spread(ti2, Genotype, nexp)
gids = ti3$gid[ti3$B73 >= 1 & ti3$Mo17 >= 1]
length(gids)

to = ti[ti$gid %in% gids,]
fo = file.path(dirw, '36.long.filtered.tsv')
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)



