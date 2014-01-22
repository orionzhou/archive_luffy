library(plyr)
library(ggplot2)

##### calculate library parameters
len_assembly = 505561755
len_bases = 376163954
len_alignable = 308662854

dir = file.path(DIR_Misc3, "hapmap_mt40/11_pipe_mapping")
f_rn = file.path(dir, "../../ncgr_fastq/04_fastq_stats.tbl")
t_rn = read.table(f_rn, header=T, sep="\t", as.is=T)[,c('sm','lb','rn','rl')]

t_lb = ddply(t_rn, .(lb), summarise, rl=unique(rl)[1])

f61 = file.path(dir, "61_stat_raw.tbl")
t61 = read.table(f61, header=T, sep="\t", as.is=T)
f62 = file.path(dir, "62_isd.tbl")
t62 = read.table(f62, header=T, sep="\t", as.is=T)


rglb = data.frame(rg=t61$rg, lb=t61$rg, stringsAsFactors=F)
for (i in 1:nrow(rglb)) {
  idxs1 = which(t_rn$rn == rglb$rg[i])
  if(length(idxs1) == 1) {
    rglb$lb[i] = t_rn$lb[idxs1[1]]
  }
}

t61.1 = merge(rglb, t61, by='rg')

ts.1 = ddply(t61.1, .(lb), summarise, 
  total=sum(total), unmapped=sum(unmapped), unpaired=sum(unpaired), unpaired_dedup=sum(unpaired_dedup), unpaired_uniq=sum(unpaired_uniq), 
  paired=sum(paired), paired_dedup=sum(paired_dedup), 
  paired_uniq=sum(paired_uniq), paired_proper=sum(paired_proper))

ts.2 = merge(t_lb, ts.1, by='lb')
ts.2$rl = as.numeric(ts.2$rl)
ts = ts.2
ts.3 = cbind(ts, 
  pct_dedup=(ts$paired_dedup * 2 + ts$unpaired_dedup) / (ts$paired * 2 + ts$unpaired), 
  cov_total=ts$total * ts$rl * 2 / len_assembly, 
  cov_dedup=((ts$paired_dedup)*2 + (ts$unpaired_dedup)) * ts$rl / len_bases,
  cov_uniq=(ts$paired_uniq * 2 + ts$unpaired_uniq) * ts$rl / len_alignable,
  is_mean=NA, is_median=NA, is_sd=NA, is_mad=NA, is_mld=NA, is_mno=NA) 


t62.1 = merge(rglb, t62, by='rg')
t62.2 = ddply(t62.1, .(lb, is), summarise, cnt=sum(cnt))

ts = ts.3
ti = t62.2
for (i in 1:nrow(ts)) {
  if(is.na(ts$total[i])) next
  lb = ts$lb[i]
  
  d1 = ti[ti$lb==lb, c('is','cnt')]
  d2 = with(d1, rep(is, cnt))
  ts$is_mean[i] = mean(d2)
  ts$is_sd[i] = sd(d2)
  ts$is_median[i] = median(d2)
  tmp = quantile(d2, c(0.159, 0.841))
  ts$is_rsd[i] = (tmp[2] - tmp[1]) / 2
  ts$is_mad[i] = mad(d2)
}
ts$is_mld = 10 * ts$is_mad
ts$is_mno = ts$is_median + 20 * ts$is_mad
write.table(ts, file.path(dir, "71_stat.tbl"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# plot coverage
s4$id = factor(s4$id, levels=sort(unique(s4$id), decreasing=TRUE))
p = ggplot(data=s4) +
  geom_bar(mapping = aes(x=id, y=rp/1000000, fill=rp_cat), stat='identity', position='stack') +
  layer(data=s3, geom="text", mapping=aes(x=acc, y=total/1000000 - 5, label=sprintf("%.01fx", cov_mapping)), geom_params=list(size=3.5, color='white')) +
  layer(data=s3, geom="text", mapping=aes(x=acc, y=total/1000000 + 5, label=sprintf("%.01fx", cov_sequencing)), geom_params=list(size=3.5)) +
  scale_fill_brewer(palette='Set1') +
  scale_x_discrete(name='') +
  scale_y_continuous(name='million read pairs') +
  coord_flip() +
  labs(fill="read pair status")
ggsave(p, filename=file.path(dirO, "cov.png"), width=8, height=6)

# plot ISD
lbs = unique(t_rg$lb[t_rg$sm %in% sms])
i3 = i2[as.character(i2$lb) %in% lbs, ]
p <- ggplot(i3, aes(is, cnt/1000000)) + 
  geom_line(aes(), size = 0.1) + 
  scale_x_continuous(name="insert size", limits=c(0,700), breaks=c(300,600)) + 
  scale_y_continuous(name="million read pairs") + 
  scale_colour_brewer(palette='Paired') +
#  opts(axis.text.x = theme_text(hjust=1, size=10)) +
  facet_wrap( ~ lb, nrow = 5)
ggsave(p, filename=file.path(dirO, "isd.png"), width=10, height=8)

#check for consistency with Joann's clonality analysis
si1 = read.table(file.path(dirW, "06_stat_isd.tbl"), sep="\t", header=TRUE)
o1 = aggregate(s2$pct_dup, by=list(s2$id), FUN=mean)
colnames(o1) = c('id', 'pct_dup')
o2 = read.table(file.path(DIR_R, "ngs/clonality.tbl"), header=TRUE, sep="\t")
o3 = merge(o1, o2, by="id", all=FALSE)
plot(o3$pct_dup, o3$pct1, type="p", xlab="proportion of duplicates", ylab="% distinct pairs occuring 1X")
fit = lm(o3$pct1 ~ o3$pct_dup)
abline(fit, col="blue")
fit.sum = summary(fit)
ann = paste("adjusted Rsquare = ", sprintf("%.04f", fit.sum$adj.r.squared), sep="")
text(0.2, 0.3, ann, col='red')

#process hydra-sv outputs
dir = file.path(DIR_Repo, "mt_35/40_sv/03_hydra")
dirO = file.path(DIR_Repo, "mt_35/40_sv/04_hydra_filtered")

for (i in 1:30) {
  id = sprintf("HM%03d", i)
  
  f01 = file.path(dir, paste(id, ".all", sep=""))
  if(file.info(f01)$size == 0) next
  f02 = file.path(dir, paste(id, ".detail", sep=""))
  d01 = read.table(f01, sep="\t", header=F, as.is=T)
  d02 = read.table(f02, sep="\t", header=F, as.is=T, comment.char="&")
  colnames(d01) = strsplit("chrom1 beg1 end1 chrom2 beg2 end2 breakpointId numDistinctPairs strand1 strand2 meanEditDist1 meanEditDist2 meanMappings1 meanMappings2 size numMappings allWeightedSupported finalSupport finalWeightedSupport numUniquePairs numAnchoredPairs numMultiplyMappedPairs blank", " ")[[1]]
  colnames(d02) = strsplit("chrom1 beg1 end1 chrom2 beg2 end2 name mate1End strand1 strand2 editDist1 editDist2 numMappings1 numMappings2 mappingType includedInBreakpoint breakpointId", " ")[[1]]
  d03 = aggregate(d02[,c('name','strand1')], by=list(factor(d02$breakpointId)), FUN=strjoin)
  colnames(d03) = c('breakpointId', 'names', 'strand1s')
  d05 = merge(d01[,1:22], d03[,c('breakpointId', 'readnames')], by=c('breakpointId'))
  
  write.table(d05, file.path(dirO, paste(id, ".tbl", sep="")), sep="\t", row.names=F, col.names=T, quote=F)
}

