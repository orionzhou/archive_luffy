require(dplyr)
require(GenomicRanges)
require(tidyr)
require(ape)
require(ggplot2)
require(plyr)

dirw = file.path(Sys.getenv("misc2"), "briggs")
dirw = '/home/springer/zhoux379/scratch/briggs'

dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
f_gtb = file.path(dirg, "51.gtb")
tg = read.table(f_gtb, sep = "\t", header = T, as.is = T)[,c('id','par','cat2','cat3')]

ft = file.path(dirg, "57.longest.tsv")
tt = read.table(ft, sep = "\t", header = F, as.is = T)
tt = tg[tg$id %in% tt$V2,]

fi = file.path(dirw, '00.1.read.correct.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)[,1:5]

### get raw read count
fi = file.path(dirw, '00.5.htseq.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)[,1:5]

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

### compute RPM
fi1 = file.path(dirw, '32.rc.tsv')
fi2 = file.path(dirw, '32.rc.as.tsv')

ti1 = read.table(fi1, header = T, sep = "\t", as.is = T)
ti2 = read.table(fi2, header = T, sep = "\t", as.is = T)

total_reads = apply(ti1[,-1]+ti2[,-1], 2, sum)

to1 = ti1
to2 = ti2
for (i in 1:length(total_reads)) {
	to1[,i+1] = to1[,i+1] / total_reads[i] * 1000000
	to2[,i+1] = to2[,i+1] / total_reads[i] * 1000000
}
fo1 = file.path(dirw, '33.rpm.tsv')
write.table(to1, fo1, sep = "\t", row.names = F, col.names = T, quote = F)
fo2 = file.path(dirw, '33.rpm.as.tsv')
write.table(to2, fo2, sep = "\t", row.names = F, col.names = T, quote = F)

### stats
fi1 = file.path(dirw, '33.rpm.tsv')
fi2 = file.path(dirw, '33.rpm.as.tsv')

ti1 = read.table(fi1, header = T, sep = "\t", as.is = T)
ti2 = read.table(fi2, header = T, sep = "\t", as.is = T)
e1 = ti1[,-1]
e2 = ti2[,-1]

tl = data.frame()
tiss = unique(ti$Tissue); gts = unique(ti$Genotype)
for (tis in tiss) {
	for (gt in gts) {
		sms = ti$SampleID[ti$Tissue == tis & ti$Genotype == gt]
		if(length(sms) == 0) next
		x1 = ti1[,sms]
		x2 = ti2[,sms]
		rpm1 = apply(x1, 1, mean)
		rpm2 = apply(x2, 1, mean)
		t1 = data.frame(gid = ti1$gid, Tissue = tis, Genotype = gt, rpm1 = rpm1, rpm2 = rpm2, stringsAsFactors = F)
		tl = rbind(tl, t1)
	}
}


#tl1 = gather(ti1, SampleID, rpm, -gid)
#tl2 = gather(ti2, SampleID, rpm, -gid)
#tl3 = rbind(cbind(tl1, type='sense'), cbind(tl2, type='antisense'))
#tl4 = merge(tl3, ti[,c(1,3:5)], by='SampleID')
#grp = dplyr::group_by(tl4, SampleID, gid, type, Tissue, Genotype)
#tl = dplyr::summarise(grp, rpm = median(rpm))

fl = file.path(dirw, '39.rpm.long.tsv')
write.table(tl, fl, sep = "\t", row.names = F, col.names = T, quote = F)


## prop of antisense transcription
sum1 = apply(ti1[,-1], 2, sum)
sum2 = apply(ti2[,-1], 2, sum)
stopifnot(names(sum1) == ti$SampleID)
tp = cbind(ti[,1:5], x = 1:nrow(ti), pct_as = sum2/(sum1+sum2))
tpx = ddply(tp, .(Tissue), summarise, x = mean(x))
p1 = ggplot(tp) +
  geom_bar(aes(x=x, y=pct_as, fill=Genotype), stat='identity') +
  scale_x_continuous(name='Tissue', breaks = tpx$x, labels = tpx$Tissue, expand = c(0,0), limits = c(0,max(tp$x)+1)) +
  scale_y_continuous(name='Proportion of AntiSense Transcription', expand = c(0.01,0)) +
  coord_flip() +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.3,1), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) + 
  geom_rect(xmin=64,xmax=65,ymin=40,ymax=41,fill='black')

fp = file.path(dirw, "stats/41.as.prop.pdf")
ggsave(p1, filename=fp, width=5, height=8)

## 
fl = file.path(dirw, '39.rpm.long.tsv')
tl = read.table(fl, header = T, sep = "\t", as.is = T)
tl$rpm1 = asinh(tl$rpm1)
tl$rpm2 = asinh(tl$rpm2)
tl1 = tl[tl$Genotype %in% c('B73', "Mo17"),]

tl21 = tl1[,c('gid','Tissue','Genotype','rpm1')] %>% spread(Genotype, rpm1)
tl22 = tl1[,c('gid','Tissue','Genotype','rpm2')] %>% spread(Genotype, rpm2)
colnames(tl22)[3:ncol(tl22)] = paste(colnames(tl22)[3:ncol(tl22)], '_as', sep = '')
tl2 = cbind(tl21, tl22[,-c(1:2)])

tl3 = tl2[(tl2$B73_as > 0 | tl2$Mo17_as > 0) & tl2$Tissue == 'sheath_v12',]


tl4 = tl3[(tl3$B73_as == 0 & tl3$Mo17_as > 1) | (tl3$B73_as > 1 & tl3$Mo17_as == 0),]
#tl5 = ddply(tl4, .(gid), summarise, ntis = length(Tissue))
#tl6 = merge(tl5, tt[,-1], by.x = 'gid', by.y = 'par')
#tl6[tl6$ntis >=2,][1:100,]

e = tl4[,-c(1:2)]
cor_opt = "pearson"
hc_opt = "ward.D"
e.r.dist <- as.dist(1-cor(t(e), method = cor_opt))
e.r.hc <- hclust(e.r.dist, method = hc_opt)
hc = e.r.hc
idxs = hc$order

tl5 = gather(tl4[,-2], type, rpm, -gid)
tl5$gid = factor(tl5$gid, levels = tl5$gid[idxs])

pb <- ggplot(tl5) +
  geom_tile(aes(x = type, y = gid, fill = rpm), height = 1) +
  theme_bw() + 
  scale_x_discrete(name = '') +
  scale_y_discrete(expand = c(0, 0), name = '') +
#  scale_fill_gradient(name = 'RPM', space = "Lab", low = 'firebrick1', high = 'dodgerblue', na.value = 'grey50') +
#  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0, 1), legend.title = element_text(size = 8), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill=NA, size=0)) +
  theme(plot.margin = unit(c(0.5,0.5,0,1), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 8, angle = 60, hjust = 1)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))
fo = sprintf("%s/stats/42.heatmap.pdf", dirw)
ggsave(pb, filename = fo, width = 8, height = 10)



### expression variability for genes w/ and w/o antisense transcripts
tl1 = tl[tl$Genotype == 'Mo17',]
grp = dplyr::group_by(tl1, gid)
tl2 = dplyr::summarise(grp, cv = sd(rpm1)/mean(rpm1)*100, ny = sum(rpm2>0))
cv1 = tl2$cv[tl2$ny>=5]
cv2 = tl2$cv[tl2$ny==0]
t.test(cv1, cv2, 'alternative' = 'less')

### check overlap of top5000 sense transcripts and anti-sense transcripts
res = c()
for (tissue in unique(tl$Tissue)) {
	for (genotype in c('B73', 'Mo17')) {
	tls = tl[tl$Tissue == tissue & tl$Genotype == genotype,]
	nsam = 1000
	gids1 = tls$gid[order(tls$rpm1)][1:nsam]
	gids2 = tls$gid[order(tls$rpm2)][1:nsam]
	res = c(res, sum(gids1 %in% gids2))
	}
}
hist(res, breaks = 10, xlim = c(0,1000))
dev.off()

res = c()
popsize = length(unique(tls$gid))
for (i in 1:1000) {
	s1 = sample(1:popsize, size = nsam)
	s2 = sample(1:popsize, size = nsam)
	res = c(res, sum(s1 %in% s2))
}
hist(res, breaks = 10, xlim = c(0,1000))
dev.off()

###
tls = tl[tl$Genotype == 'B73',]
grp = dplyr::group_by(tls, gid)
t1 = dplyr::summarise(grp, b_n2 = sum(rpm2>0), b_pcc = cor(rpm1, rpm2), b_pval = cor.test(rpm1, rpm2)$p.value)

tls = tl[tl$Genotype == 'Mo17',]
grp = dplyr::group_by(tls, gid)
t2 = dplyr::summarise(grp, m_n2 = sum(rpm2>0), m_pcc = cor(rpm1, rpm2), m_pval = cor.test(rpm1, rpm2)$p.value)
t2 = merge(t1, t2, by = 'gid')

tls = tl[tl$Genotype == 'B73xMo17',]
grp = dplyr::group_by(tls, gid)
t3 = dplyr::summarise(grp, h_n2 = sum(rpm2>0), h_pcc = cor(rpm1, rpm2), h_pval = cor.test(rpm1, rpm2)$p.value)
t3 = merge(t2, t3, by = 'gid')

t4 = t3[which((t3$b_n2 >= 10 & t3$b_pcc < 0 & t3$b_pval < 0.05) | (t3$m_n2 >= 10 & t3$m_pcc < 0 & t3$m_pval < 0.05)),]
nrow(t4) # 242
t41 = t4[t4$b_pval < 0.05 & t4$m_pval < 0.05,]
nrow(t41) # 81
sum(t41$h_pval < 0.05, na.rm = T) # 54  67%
t42 = t4[which(!(t4$b_pval < 0.05 & t4$m_pval < 0.05)),]
nrow(t42) # 161
sum(t42$h_pcc < 0 & t42$h_pval < 0.05) # 53 33%


t8 = t3[which((t3$b_n2 >= 10 & t3$b_pcc > 0 & t3$b_pval < 0.05) | (t3$m_n2 >= 10 & t3$m_pcc > 0 & t3$m_pval < 0.05)),] # 7697
t81 = t8[t8$b_pval < 0.05 & t8$m_pval < 0.05,] #4424
sum(t81$h_pval < 0.05, na.rm = T) # 3738  84%
t82 = t8[which(!(t8$b_pval < 0.05 & t8$m_pval < 0.05)),] # 3273
sum(t82$h_pcc > 0 & t82$h_pval < 0.05, na.rm = T) # 1491 45%



