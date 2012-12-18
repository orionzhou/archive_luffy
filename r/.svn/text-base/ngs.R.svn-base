sms = get_mt_ids("acc33")
len_assembly = 505561755
len_bases = 376163954
len_alignable = 308662854

dir = "/project/youngn/zhoup/Data/misc3/hapmap"
dirO = "/project/youngn/zhoup/Data/misc3/hapmap/80_stat"
f_rg = file.path(dir, "09_fastq.tbl")
t_rg = read.table(f_rg, header=T, sep="\t", as.is=T)[,c('sm','rg','lb','rl','pi')]

f13a = file.path(dir, "15_pipe_bam/13_stat.tbl")
s1 = read.table(f13a, header=T, sep="\t", as.is=T)
s1.1 = merge(t_rg[,c('rg','lb')], s1, by='rg')
s1.2 = s1.1[,-1]

sum_by_col <- function(df, cols) {
  if(length(cols) == 1) {
    sum(df[,cols])
  } else {
    colSums(df[,cols])
  }
}
s1.3 = ddply(s1.2, .(lb), sum_by_col, cols=2:10)
s1.4 = merge(unique(t_rg[,c('lb','rl')]), s1.3)

s2 = cbind(s1.4, pct_dup=NA, cov_total=NA, cov_rmdup=NA, cov_uniq=NA, is_mean=NA, is_median=NA, is_sd=NA, is_mad=NA, is_mld=NA, is_mno=NA)
s2$pct_dup = (s2$paired_dup * 2 + s2$unpaired_dup) / (s2$paired * 2 + s2$unpaired)
s2$cov_total = s2$total * s2$rl * 2 / len_assembly 
s2$cov_rmdup = ((s2$paired-s2$paired_dup)*2 + (s2$unpaired-s2$unpaired_dup)) * s2$rl / len_bases 
s2$cov_uniq = (s2$paired_uniq * 2 + s2$unpaired_uniq) * s2$rl / len_alignable


f13b = file.path(dir, "15_pipe_bam/13_stat_isd.tbl")
i1 = read.table(f13b, header=T, sep="\t", as.is=T)
i1.1 = merge(t_rg[,c('rg','lb')], i1, by='rg')
i1.2 = i1.1[,-1]

i2 = ddply(i1.2, .(lb, is), sum_by_col, cols=3)
colnames(i2)[3] = "cnt"

for (i in 1:nrow(s2)) {
  if(is.na(s2$total[i])) next
  d1 = i2[i2$lb == s2$lb[i], c('is','cnt')]
  d2 = with(d1, rep(is, cnt))
  s2$is_mean[i] = mean(d2)
  s2$is_sd[i] = sd(d2)
  s2$is_median[i] = median(d2)
  tmp = quantile(d2, c(0.159, 0.841))
  s2$is_rsd[i] = (tmp[2] - tmp[1]) / 2
  s2$is_mad[i] = mad(d2)
}
s2$is_mld = 10 * s2$is_mad
s2$is_mno = s2$is_median + 20 * s2$is_mad
write.table(s2, file.path(dir, "15_pipe_bam/14_stat.tbl"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

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

