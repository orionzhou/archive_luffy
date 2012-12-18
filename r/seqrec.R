dir = file.path(DIR_Misc1, "seq07")
f_info = file.path(dir, "32_info.tbl")
d01 = read.table(f_info, header=T, sep="\t", stringsAsFactors=T)
d02 = d01[d01$type == "cds",]

gf1 = read.table(file.path(DIR_Misc2, "genefam/41_genefam.tbl"), sep="\t", header=T, as.is=T, quote="")

pre = "21_seq"
accs = get_mt_ids("acc26")
d03 = d02[d02$acc %in% accs,]
d04 = cbind(d03[,c(-3)], cov=1-d03$num_N/d03$length, snpDensity=d03$num_snp/(d03$length-d03$num_N))
d05 = merge(d04, gf1[,c(1:3)], by="id")
p = ggplot(d05, aes(x=type2, y=cov, fill=type2)) +
  geom_boxplot() + 
  scale_fill_brewer(palette='Set1') +	
  labs(fill="Gene Family") +
  scale_y_continuous(name='% coverage') +
  scale_x_discrete(name='') +
  facet_wrap( ~ acc, ncol=9) +
  opts(axis.text.x=theme_blank())
ggsave(p, filename = file.path(dir, pre, "31_cov.png"), width=12, height=9)
p = ggplot(d05, aes(x=type2, y=snpDensity, fill=type2)) +
  geom_boxplot() + 
  scale_fill_brewer(palette='Set1') +	
  labs(fill="Gene Family") +
  scale_y_continuous(name='SNP density', limits=c(0, 0.06)) +
  scale_x_discrete(name='') +
  facet_wrap( ~ acc, ncol=9) +
  opts(axis.text.x=theme_blank())
ggsave(p, filename = file.path(dir, pre, "32_snpDensity.png"), width=12, height=9)


pre = "22_seq"
accs = get_mt_ids("opt12")
outgroups = c("HM102")
napi_thresh = 40
napo_thresh = 1
nsm_thresh = 1

pre = "23_seq"
accs = get_mt_ids("opt13")
outgroups = c("HM330")
napi_thresh = 40 
napo_thresh = 1
nsm_thresh = 1

d03 = d02[d02$acc %in% accs,]
d04 = cbind(d03, cov=1-d03$num_N/d03$length)
d11i = d04[!d04$acc %in% outgroups,]
d11o = d04[d04$acc %in% outgroups,]

cov_thresh = 0.7
cnt_acc_passed <- function(a) sum(a>=cov_thresh)
r01 = aggregate(d04$length, by=list(factor(d04$id)), FUN=unique)
colnames(r01) = c("id", "length")
r11 = aggregate(d11i$cov, by=list(factor(d11i$id)), FUN=cnt_acc_passed)
r12 = aggregate(d11o$cov, by=list(factor(d11o$id)), FUN=cnt_acc_passed)
r21 = aggregate(d11i$num_snp, by=list(factor(d11i$id)), FUN=max)
r22 = aggregate(d11o$num_snp, by=list(factor(d11o$id)), FUN=max)

r31 = cbind(r01[order(r01$id),], num_acc_passed_ingroup=r11[order(r11[,1]), 2], num_acc_passed_outgroup=r12[order(r12[,1]), 2], num_snp_max_ingroup=r21[order(r21[,1]), 2], num_snp_max_outgroup=r22[order(r22[,1]), 2])
r32 = cbind(r31, num_acc_passed=r31$num_acc_passed_ingroup+r31$num_acc_passed_outgroup, num_snp_max=apply(r31[,c(5,6)], 1, max))
r33 = cbind(r32, select=as.numeric(r32$num_acc_passed_ingroup >= napi_thresh & r32$num_acc_passed_outgroup>=napo_thresh & r32$num_snp_max>=nsm_thresh))

f01 = file.path(dir, pre, "01_stat.tbl")
write.table(r33, f01, sep="\t", col.names=T, row.names=F, quote=F)
f02 = file.path(dir, pre, "02_ids.tbl")
write.table(r33[r33$select==1, c('id','length')], f02, sep="\t", col.names=T, row.names=F, quote=F)


