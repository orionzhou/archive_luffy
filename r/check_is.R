library(lattice)
fi = '/home/youngn/zhoup/Data/misc3/ncgr_fastq/04_fastq_stats.tbl'
fq = read.table(fi, header=T, as.is=T, quote='')

d11 = '/home/youngn/zhoup/Data/misc3/hapmap_mt40/11_pipe_mapping'
sm = "HM056"
lbs = unique(fq$lb[fq$sm==sm])

for (i in 1:length(lbs)) {
	lb = lbs[i]
	fi = sprintf("%s/12_dedup/%s.isd.tbl", d11, lb)
	ti = read.table(fi, header=T, sep="\t", as.is=T)
	ti = cbind(lb=lb, ti)
	if(i == 1) {
		t_is = ti
	} else {
		t_is = rbind(t_is, ti)
	}
}
xyplot( cnt ~ is | lb, data = t_is, type='l', auto.key=T)

lb = "HM340_LIPE_1"
lb = "HM056_LIPE_1"
pres = c("12_dedup", "19_remap_dedup")
for (i in 1:length(pres)) {
	pre = pres[i]
	fi = sprintf("%s/%s/%s.isd.tbl", d11, pre, lb)
	ti = read.table(fi, header=T, sep="\t", as.is=T)
	ti = cbind(pre=pre, ti)
	if(i == 1) {
		t_is = ti
	} else {
		t_is = rbind(t_is, ti)
	}
}
xyplot( cnt ~ is, data = t_is, type='l', ylim=c(0,1000), groups=pre, auto.key=T)
