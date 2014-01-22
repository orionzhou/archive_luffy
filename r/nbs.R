org = "Mtruncatula_4.0"
dir = sprintf("%s/spada.nbs.%s/31_model_evaluation", DIR_Misc4, org)
fs = file.path(dir, "41_stat.tbl")
fg = file.path(dir, "55_nonovlp.gtb")
ts = read.table(fs, head=T, sep="\t")[,1:14]
tg = read.table(fg, head=T, sep="\t")[,1:6]

tg2=merge(tg,ts,by="id")
tgs = tg2[tg2$e<=0.001 & tg2$score_aln>=-1000,]
tgs = tgs[order(tgs$e),]
nrow(tgs)
