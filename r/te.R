source("Location.R")

dirg = file.path(Sys.getenv("genome"), "HM101")

fr = file.path(dirg, "12.rm.tbl")
tr = read.table(fr, sep = "\t", header = F, as.is = T)
#idxs = grep("^(DNA)|(LINE)|(LTR)|(SINE)|(RC)", tr$V9)
grr = with(tr, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

fi = file.path(dirg, "raw/41.gtb")
ti = read.table(fi, sep = "\t", header = T, as.is = T, quote = "")[,c(1:6,16:18)]
gri = with(ti, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

lens = intersect_basepair(gri, grr)
idxs_rm_pfam = lens/(ti$end-ti$beg+1) >= 0.6 | ti$cat2 == 'TE'
idxs_rm_pfam = ti$cat2 == 'TE'
idxs = grep("Medtr[0-9]+te", ti$id)
idxs_jcvi = rep(FALSE, nrow(ti))
idxs_jcvi[idxs] = TRUE

table(data.frame(jcvi=idxs_jcvi,rm_pfam=idxs_rm_pfam))


fj = file.path(dirg, "raw/25.dedup.gtb")
tj = read.table(fj, sep = "\t", header = T, as.is = T, quote = "")[,c(1:6,16:18)]

tis = ti[ti$cat2 %in% c("Metallophos", "PTR2", "DIOX_N", "B3", "RRM_1"),]

idxs = grep("Medtr[0-9]+te", tis$id)

x = table(tis$cat2)
names(x)[1] = 'none'
y = table(tis$cat2[idxs])
names(y)[1] = 'none'

tf = data.frame(fam = names(x), cnt = as.numeric(x), stringsAsFactors = F)
tf = cbind(tf, cntt = y[tf$fam])
tf = cbind(tf, pct = tf$cntt / tf$cnt)
tf = tf[order(tf$pct, tf$cntt, decreasing = T),]
tf[1:100,]



### check ovlp of Mt4.0 TEs w. RepeatMasker
t1 = ti[idxs,]
t1 = t1[t1$cat2 == 'Unknown',]
gr1 = with(t1, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

bps = intersect_basepair(gr1, grr)
pct = bps / (t1$end - t1$beg + 1)
hist(pct)

DUF4283 bZIP_1 DUF4371
tj[tj$id %in% ti$id[ti$cat3 == 'DUF4371'], c('id','note')]

