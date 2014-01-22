#read coverage statistics
f_cov_stat = file.path(DIR_Repo, "mt_35/30_vnt_acc84/21_coverage/11_info.tbl")
covStat = read.table(f_cov_stat, header=T, sep="\t")
acc1 = "HM018"
acc2 = "HM101"
cov2_ratio_adjust = covStat$cov2_sum[covStat$acc==acc1 & covStat$chr=="chr5"] / covStat$cov2_sum[covStat$acc==acc2 & covStat$chr=="chr5"]

#compare aCGH ratio with reseq-ratio for SV regions
f_cov = file.path(DIR_Misc2, "acgh/52_cov.tbl");
cov = read.table(f_cov, header=T, sep="\t")

dirW = file.path(DIR_Misc2, "acgh", acc1)
cov01 = cov[cov$acc==acc1,]
cov02 = cov[cov$acc==acc2,]
cov03 = merge(cov01, cov02, by=c("id", "chr", "beg", "end"))
cov04 = cbind(cov03, cov2_ratio = log2((cov03$cov2_mean.x/cov03$cov2_mean.y)/cov2_ratio_adjust))
write.table(cov04, file.path(dirW, "09_stat.tbl"), sep="\t", row.names=F, col.names=T, quote=F)

f_sv = file.path(dirW, "23_sv_gene.tbl")
g01 = read.table(f_sv, header=T, sep="\t", stringsAsFactors=F)
ids = unlist(strsplit(g01$idG, " "))

cov_all = cov04$cov2_ratio
cov_all = cov_all[!is.na(cov_all) & is.finite(cov_all)]
cov_s = cov04$cov2_ratio[cov04$id %in% ids]
cov_s = cov_s[!is.na(cov_s) & is.finite(cov_s)]
cov_us = cov04$cov2_ratio[!cov04$id %in% ids]
cov_us = cov_us[!is.na(cov_us) & is.finite(cov_us)]
t.test(cov_s, cov_us, alternative="less")

n_all = length(cov_all)
n_called = sum(cov_s < quantile(cov_all, 0.05))

d01 = data.frame(cov_s = cov_s)
d02 = data.frame(cov_us = cov_us)
p <- ggplot(data=d02) +
  geom_histogram(mapping=aes(x=cov_us, fill='genes NOT in SV regions'), geom_params=list(alpha=0.5)) +
  layer(data=d01, geom='histogram', mapping=aes(x=cov_s, fill='genes IN SV regions'), geom_params=list(alpha=1)) +
  scale_fill_brewer(palette="Set1", name='') +
  scale_x_continuous(name='log2(average_unique_coverage_ratio)') +
  scale_y_continuous(name='gene count') 
ggsave(p, filename = sprintf("%s/52_%s.png", dirW, acc1), width=5, height=7)

#check for correlation of nimblegen arrayCGH results and resequencing depth
dirW = file.path(DIR_Misc2, "acgh")
fseq = file.path(dirW, "08_cov.tbl")
s1 = read.table(fseq, header=TRUE, sep="\t")

fcgh = file.path(dirW, acc1, "10.tbl")
c1 = read.table(fcgh, header=TRUE, sep="\t")
c2 = c1[,c("PROBE_ID", "RATIO_CORRECTED")]
colnames(c2) = c("id", "ratio_acgh")

f_probe_1 = file.path(dirW, acc1, "21_sv_probes.tbl")
probes = read.table(f_probe_1, header=T, sep="\t")$Probe.ID

s1.1 = s1[s1$acc == acc1, ] 
s1.2 = s1[s1$acc == acc2, ] 
s2 = merge(s1.1, s1.2, by=c("id", "chr", "beg", "end"))

s3 = s2[!s2$id %in% probes, ]
s4 = cbind(s3, cov2_ratio = log2((s3$cov2_mean.x/s3$cov2_mean.y)/cov2_ratio_adjust))
s5 = s4[s4$cov2_mean.x>5 & s4$cov2_sd.x/s4$cov2_mean.x<0.1 & s4$cov2_mean.y>5 & s4$cov2_sd.y/s4$cov2_mean.y<0.1, ]

d1 = merge(c2, s5, by="id", all=FALSE)
d2 = d1[is.finite(d1$cov2_ratio),]
f_png1 = file.path(dirW, acc1, "51.png")
lg(d2$ratio_acgh, d2$cov2_ratio, "log2-ratio(aCGH)", "log2-coverage-ratio(resequencing)", f_png1)



