dirW = "/home/youngn/zhoup/Data/genome/Mtruncatula_4.0"
fc = file.path(dirW, "51_seqlen.tbl")
tc = read.table(fc, sep="\t", header=TRUE, as.is=TRUE)
tc = tc[tc$id < 'chr9',]

fg = file.path(dirW, "52_gap_loc.tbl")
tg = read.table(fg, sep="\t", header=TRUE, as.is=TRUE)
tg = tg[tg$id < 'chr9',]

tt = table(tg$len)
x1 = data.frame(gap_len=as.numeric(names(tt)), gap_cnt=as.numeric(tt))
x2 = cbind(x1, gap_cumsum=cumsum(x1$gap_len*x1$gap_cnt))
p <- ggplot(data=x2) +
  geom_line(mapping=aes(x=gap_len, y=gap_cumsum/1000000), geom_params=list(col='blue')) +
  scale_x_log10(name='Gap size (bp)', breaks = 10^(1:5), labels=sprintf("10^%d", (1:5))) +
  scale_y_continuous(name='Cumulative sum of gap contribution (Mbp)') 
ggsave(p, filename = file.path(dirW, "gap_stat.png"), width=6, height=4)

tg = tg[tg$len >= 1000,]
tt = tapply(tg$len, as.factor(tg$id), sum)
ttl = data.frame(id=names(tt), total_gap_len=tt)
tt = tapply(tg$len, as.factor(tg$id), length)
ttn = data.frame(id=names(tt), total_gap_count=tt)
tt = tapply(tg$len, as.factor(tg$id), mean)
ttm = data.frame(id=names(tt), mean_gap_len=tt)
tt = tapply(tg$len, as.factor(tg$id), median)
ttd = data.frame(id=names(tt), median_gap_len=tt)
ts = merge(tc[,c('id','length')], ttl, by='id')
ts = merge(ts, ttn, by='id')
ts = merge(ts, ttm, by='id')
ts = merge(ts, ttd, by='id')

write.table(ts, file=file.path(dirW, "52_gap_stat.tbl"), col.names=T, row.names=F, sep="\t", quote=F)

ts = t1[t1$id > 'chr9' & t1$type == 'chr',]
