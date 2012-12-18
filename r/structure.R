dirW = file.path(DIR_Misc3, "structure/12_acc261")
rep = 1

  f01 = file.path(dirW, "03_prob.txt");
  stat = read.table(f01, header=TRUE);
  p <- ggplot(stat, aes(x=factor(stat$K), y=Mean.P.D.., ymin=Mean.P.D..-Var.P.D.., ymax=Mean.P.D..+Var.P.D.., color=factor(stat$Run))) + 
    geom_pointrange(position=position_dodge(width=0.5)) +
#    scale_colour_brewer(palette="Set1") +
    scale_x_discrete(name='K') +
    scale_y_continuous(name='P(D|K)') +
    labs(colour="Run") +
    opts(title='');
  ggsave(p, filename = file.path(dirW, "03_prob.png"), width=7, height=5);

for (k in 2:10) {
  fi = sprintf("%s/05_stat/k_%d_rep_%d.tbl", dirW, k, rep)
  d01 = read.table(fi, header=T, sep="\t", stringsAsFactors=F)

  fp = sprintf("%s/06_figs/k_%d_rep_%d.png", dirW, k, rep)
  p = ggplot(d01) +
    geom_bar(aes(x=id, y=mean, fill=factor(k)), stat='identity', position='stack') +
    scale_fill_brewer(palette='Set3', name='ancestry') +  
    scale_x_discrete(name='', breaks=sprintf("HM%03d", c(1,60,101,317))) + 
    scale_y_continuous(name='') +
    opts(axis.text.x = theme_text(hjust=0, size=10)) +
    opts(axis.text.y = theme_blank())
  ggsave(p, filename = fp, width=10, height=4)
}

