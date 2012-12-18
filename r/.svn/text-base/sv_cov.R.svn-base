#read gene annotation
f_ann = file.path(DIR_Genome, "mt_35/10_model_Mt3.5v6/86.tbl")
a01 = read.table(f_ann, header=T, as.is=T, sep="\t", quote="")

#read coverage statistics
f_cov_stat = file.path(DIR_Repo, "mt_35/30_vnt_acc84/21_coverage/11_info.tbl")
covStat = read.table(f_cov_stat, header=T, sep="\t")
covStat = cbind(covStat, cov1_mean = covStat$cov1_sum/covStat$pos_sum, cov2_mean = covStat$cov2_sum/covStat$pos_sum)

accs = get_mt_ids("acc26")
dir = file.path(DIR_Repo, "mt_35/40_sv/41_shared")
s1 = read.table(file.path(dir, "11.tbl"), sep="\t", as.is=T, header=T, quote="", comment.char="/")
g1 = read.table(file.path(dir, "71_genotype.tbl"), sep="\t", as.is=T, header=T, quote="", comment.char="/")

d_bam = file.path(DIR_Repo, "mt_35/15_pipe_bam/01_pos_sorted")
f_bams = paste(d_bam, "/", accs, ".bam", sep="")
names(f_bams) = accs

alleleMapping = c("A"="Ref", "G"="Alt", "N"="N")

dirO = file.path(DIR_Repo, "mt_35/40_sv/90_stat")
for (i in 1:nrow(s1)) {
  plot_sv_cov(i, s1, f_bams, accs, alleleMapping, covStat, dirO)
  if(s1$end[i] - s1$beg[i] + 1 > 1000) next
  cat(i, s1$id_pindel[i], 'done\n', sep=' ')
}

plot_sv_cov <- function(i, s1, f_bams, accs, alleleMapping, covStat, dirO) {
  sv_id = s1$id_pindel[i]
  chr = s1$chr[i]
  sv_beg = s1$beg[i]
  sv_end = s1$end[i]
  sv_len = sv_end - sv_beg + 1
  beg = sv_beg - sv_len
  end = sv_end + sv_len

  cov01 = data.frame(matrix(NA, ncol=3, nrow=0))
  colnames(cov01) = c('acc', 'pos', 'coverage')
  for (acc in accs) {
    cov_raw = get_bam_coverage(f_bams[acc], chr, beg, end)
    cov_normalized = cov_raw / covStat$cov2_mean[covStat$acc==acc & covStat$chr==chr]
    cov01 = rbind(cov01, data.frame(acc=acc, pos=1:(end-beg+1), coverage=cov_normalized))
  }
  geno = data.frame(acc=colnames(g1)[-1], allele=alleleMapping[as.character(g1[g1$id==sv_id,-1])])
  cov02 = merge(cov01, geno, by='acc')
  cov02$allele = factor(cov02$allele, levels=c("Gene", "Ref", "Alt", "N"))

  a11 = a01[a01$chr == chr & a01$beg <= end & a01$end >= beg,]
  a12 = data.frame(chr=a11$chr, beg=a11$beg-beg+1, end=a11$end-beg+1, type=a11$type, id=a11$id, allele="Gene")

  ymax = max(cov02$coverage)
  ybreaks = as.numeric(sprintf("%.01f", seq(ymax/4, by=ymax/4, length=3)))
  p <- ggplot(cov02) +
    geom_rect(aes(xmin=sv_len, xmax=sv_len*2, ymin=0, ymax=max(cov02$coverage)), fill='snow2') + 
    geom_line(aes(x=pos, y=coverage, color=acc), size=0.3) +
    layer(data=a12, geom="rect", mapping=aes(xmin=beg-1, xmax=end, ymin=0, ymax=max(cov02$coverage)/9, fill=type), geom_params=list(size=0)) +
    facet_grid(allele ~ . , scale='free', space='free') +
    scale_colour_manual(values=rainbow(26), legend=FALSE) +
    scale_fill_manual(values=c('cds'='royalblue', 'intron'='oldlace', 'utr5'='burlywood', 'utr3'='tan'), legend=TRUE) +
    scale_x_continuous(name="", formatter='comma', breaks=seq(0,by=sv_len,length=4)) +
    coord_cartesian(xlim = c(0, end-beg), wise=TRUE) +
    scale_y_continuous(name="Normalized Coverage", breaks=ybreaks) +
    theme_bw() +
    opts(title=sprintf("Deletion (id=%s, location=[%d-%d], length=%dbp)\n", sv_id, sv_len+1, 2*sv_len, sv_len)) +
    opts(strip.text.y=theme_text(angle=0)) +
    labs(fill='Gene structure') +
#    opts(legend.title=theme_blank(), legend.text=theme_blank(), legend.key=theme_blank()) + 
    opts(axis.title.y=theme_text(angle=85))
  ggsave(file.path(dirO, "11_cov", paste(sv_id, ".png", sep="")), p, width=6, height=4)
}




