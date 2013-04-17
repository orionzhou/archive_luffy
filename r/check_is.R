library(lattice)
d11 = '/home/youngn/zhoup/Data/misc3/hapmap_mt40/11_pipe_mapping'
lbs = data.frame(
	lb = c("hm340.lipe.1.1", "hm340.lipe.1.3", "hm340.lipe.1.5"),
	id = c("hm340.original", "hm340.a10000", "hm340.remap")
)

lbs = data.frame(
	lb = c("hm056.lipe.1.1", "hm056.lipe.1.5"),
	id = c("hm056.original", "hm056.remap")
)

for (i in 1:nrow(lbs)) {
	fi = sprintf("%s/%s.isd.tbl", d11, lbs$lb[i])
	ti = read.table(fi, header=T, sep="\t", as.is=T)
	ti = cbind(id=lbs$id[i], ti)
	if(i == 1) {
		t_is = ti
	} else {
		t_is = rbind(t_is, ti)
	}
}
xyplot( cnt ~ is, data = t_is, type='l', groups=id, auto.key=T)