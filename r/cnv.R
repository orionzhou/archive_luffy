a1 = readAssembly("mt_35")
a = drop.levels(subset(a1[a1$chr == "chr5",], type!='centromere'))
c = drop.levels(subset(a1[a1$chr == "chr5",], type=='centromere'))

dir = file.path(DIR_Misc2, "cnv")
acc1 = "HM010"
acc2 = "HM004"
accR = "HM101"
chr = "chr5"

f_cov_stat = file.path(DIR_Repo, "mt_35/30_vnt_acc288/21_coverage/11_info.tbl")
covStat = read.table(f_cov_stat, header=T, sep="\t")
nf1 = covStat$cov2[covStat$acc==acc1 & covStat$chr==chr]
nf2 = covStat$cov2[covStat$acc==acc2 & covStat$chr==chr]
nfR = covStat$cov2[covStat$acc==accR & covStat$chr==chr]

dirI = file.path(DIR_Repo, "mt_35/30_vnt_acc288/21_coverage/02/coverage_unique")
v1 = scan(file.path(dirI, acc1, chr), what=integer(0))
v2 = scan(file.path(dirI, acc2, chr), what=integer(0))
vR = scan(file.path(dirI, accR, chr), what=integer(0))

f_win = file.path(dir, "01_windows.tbl")
winSize = 2000
w1 = read.table(f_win, sep="\t", stringsAsFactors=F, header=T)
w2 = w1[w1$end-w1$beg+1 == winSize,]

subWinSize = 100
n_subwindow = winSize/subWinSize
d1 = cbind(w2, mean1=NA, mean2=NA, meanR=NA, median1=NA, median2=NA, medianR=NA, p1R=NA, p2R=NA, p12=NA)
for(i in 1:nrow(d1)) {
  beg = d1$beg[i]
  end = d1$end[i]
  covs1 = v1[beg:end]
  covs2 = v2[beg:end]
  covsR = vR[beg:end]
  d1[i, c("mean1", "mean2", "meanR", "median1", "median2", "medianR")] = c(mean(covs1), mean(covs2), mean(covsR), median(covs1), median(covs2), median(covsR))

  sample1 = covs1[50+(1:n_subwindow-1)*100] / nf1
  sample2 = covs2[50+(1:n_subwindow-1)*100] / nf2 
  sampleR = covsR[50+(1:n_subwindow-1)*100] / nfR
  d1$p1R[i] = t.test(sample1, sampleR, paired=T, alternative="two.sided")$p.value
  d1$p2R[i] = t.test(sample2, sampleR, paired=T, alternative="two.sided")$p.value
  d1$p12[i] = t.test(sample1, sample2, paired=T, alternative="two.sided")$p.value
} 

d2 = d1[d1$median1 >= 5 | d1$median2 >=5 | d1$medianR >=5,]
d3 = cbind(d2, p1R_j=p.adjust(d2$p1R, method="bonferroni"), p2R_j=p.adjust(d2$p2R, method="bonferroni"), p12_j=p.adjust(d2$p12, method="bonferroni"))
d4 = cbind(d3, sig1R=as.numeric(d3$p1R_j<0.05), sig2R=as.numeric(d3$p2R_j<0.05), sig12=as.numeric(d3$p12_j<0.05))
d5 = cbind(d4, c1R = d4$sig1R * (as.numeric(d4$mean1/nf1>d4$meanR/nfR)+1), c2R = d4$sig2R * (as.numeric(d4$mean2/nf2>d4$meanR/nfR)+1), c12 = d4$sig12 * (as.numeric(d4$mean1/nf1>d4$mean2/nf2)+1))
write.table(d5, file.path(dir, "31.tbl"), sep="\t", row.names=F, col.names=T, quote=F)

d5 = read.table(file.path(dir, "31.tbl"), sep="\t", head=T, stringsAsFactors=F)

d6.1 = cbind(d5[,c('beg', 'end')], p=-log10(d5$p1R_j), col=d5$c1R, covr=log2((d5$mean1/d5$meanR)/(nf1/nfR)))
d6.2 = cbind(d5[,c('beg', 'end')], p=-log10(d5$p2R_j), col=d5$c2R, covr=log2((d5$mean2/d5$meanR)/(nf2/nfR)))
d6.3 = cbind(d5[,c('beg', 'end')], p=-log10(d5$p12_j), col=d5$c12, covr=log2((d5$mean1/d5$mean2)/(nf1/nf2)))
d = d6.2
fp="12.png"
fig_title = paste(acc2, ":", accR, sep="")
p = ggplot(data = d) + 
  geom_point(mapping=aes(x=(beg+end)/2, y=covr, color=factor(col)), size=1) +
  layer(data=a, geom='rect', mapping=aes(xmin=start, xmax=end, ymin=-0.2, ymax=0, fill=type), geom_params=list(size=0)) +
  layer(data=c, geom='point', mapping=aes(x=(start+end)/2, y=-0.5, shape=type), geom_params=list(size=2)) +
  scale_color_manual(name=fig_title, values=c('0'='lightcyan1', '1'='red', '2'='blue'), breaks=c(0,1,2), labels=c("insignificant", "loss", "gain")) +
  scale_fill_manual(values=c('phase 1'='aquamarine', 'phase 3'='aquamarine4'), legend=TRUE) +
  scale_shape_manual(values=c('centromere'=17), legend=TRUE) +
  labs(colour="Density", fill="Assembly", shape="") +
  scale_x_continuous(name='chr position (bp)', formatter="comma") +
  scale_y_continuous(name='log2(coverage ratio)') +
  opts(axis.text.y=theme_text(hjust=1, size=8, colour="blue"))
ggsave(p, filename = file.path(dir, fp), width=12, height=4)

str = paste(as.numeric(d5$c12==1), sep="", collapse="")
rst = gregexpr("1+", str)
idxs = data.frame(idx=rst[[1]], length=attr(rst[[1]], 'match.length'))
table(idxs$length)


