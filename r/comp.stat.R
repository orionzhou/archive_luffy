source("comp.plot.fun.R")

#acc = "hm056"
acc = "hm340"
dir = file.path(DIR_Misc3, acc)

f_len = file.path(dir, "11_seqlen.tbl")
t_len = read.table(f_len, header=TRUE, sep="\t", as.is=T)

f_gap = file.path(dir, "12_gaploc.tbl")
t_gap = read.table(f_gap, header=TRUE, sep="\t", as.is=T)

f_aln = file.path(dir, "21_blastn/05_tiled.tbl")
t_aln = read.table(f_aln, header=TRUE, sep="\t", as.is=T)

#####calculate assembly statistics
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/contigStats.R") 
library(Biostrings)
assembly <- readDNAStringSet(sprintf("%s/01_assembly.fa", dir, acc), "fasta")
N <- list(acc=width(assembly))
reflength <- sapply(N, sum)
stats <- contigStats(N=N, reflength=reflength, style="data")
stats[["Contig_Stats"]]

##### plot scaffold size distribution
tmp = cut(t_len$length/1000, breaks=c(0,1,5,10,50,100,500,1000,5000))
p = ggplot(data.frame(size=tmp)) +
  geom_bar(aes(x=factor(size)), width=0.7) + 
  scale_x_discrete(name="Scaffold Size (kb)") + 
  scale_y_continuous(name="") +
  theme(axis.text.x = element_text(angle=15, size=8))
ggsave(file.path(dir, "figs/01_scaffold_size.png"), p, width=5, height=4)

##### plot genome distribution
p <- ggplot(data=tw) +
  geom_rect(mapping=aes(xmin=hBeg/1000000, xmax=hEnd/1000000, ymin=0, ymax=1, fill=hId)) +  
  layer(data=tg, geom='rect', mapping=aes(xmin=hBeg/1000000, xmax=hEnd/1000000, ymin=-1, ymax=0), geom_params=list()) +
  scale_x_continuous(name='Chr Position (Mbp)', expand=c(0.01, 0)) +
  scale_y_continuous(name='', expand=c(0.04, 0)) +
  facet_grid(hId ~ .) + 
  theme(legend.position='right', legend.title=element_blank()) +
  theme(axis.text.x=element_text(size=8, angle=0)) +
  theme(axis.text.y=element_blank(), axis.ticks=element_blank())
ggsave(p, filename=file.path(dir, "figs/03_coverage.png"), width=7, height=5)

##### assign alignment blocks
list.tmp = assign_block_mapping(t_aln)
t11 = list.tmp$df
t12 = list.tmp$dfb
t12b = cbind(t12, qGap=t12$qLen-t12$qLen_aln, hGap=t12$hLen-t12$hLen_aln)
t12c = t12b[t12b$hGap > 100000,]

t13 = merge(t12, t_len, by='qId')
t14 = cbind(t13, pct_cov=t13$qLen_aln/t13$len_scaf)

sum_tw <- function(cutoff, df) {
  df = df[df$pct_cov >= cutoff,]
  c('num_qId'=length(unique(df$qId)), 'num_block'=nrow(df), 'qLen_aln'=sum(df$qLen_aln), 'hLen_aln'=sum(df$hLen_aln), 'qLen'=sum(df$qLen))
}
ldply(seq(0,0.3,0.05), sum_tw, t14)

t15 = t14[t14$pct_cov >= 0.05,]

f_aln = file.path(dir, '21_blastn', "15_aln.tbl")
write.table(t11, file=f_aln, col.names=T, row.names=F, sep="\t", quote=F)

f_blk = file.path(dir, '21_blastn', "16_blk.tbl")
write.table(t15, file=f_blk, col.names=T, row.names=F, sep="\t", quote=F)

##### read in alignment/block table
f_blk = file.path(dir, '21_blastn', "16_blk.tbl")
t_blk = read.table(f_blk, header=TRUE, sep="\t", as.is=T)

f_aln = file.path(dir, '21_blastn', "15_aln.tbl")
t_aln = read.table(f_aln, header=TRUE, sep="\t", as.is=T)

##### get unmapped scaffolds for NR blast
f05 = file.path(dir, "21_blastn/05_tiled.tbl")
t05 = read.table(f05, header=TRUE, sep="\t", as.is=T)
t05  = t05[order(t05$qId, t05$qBeg),]

t21 = ddply(t05, .(qId), summarise, qLen=sum(qLen), hLen=sum(hLen), n_seg=length(hId), pct=mean(pct), e=min(e), score = sum(score))
t21b = merge(t21, t_len, by="qId", all=T)
t21b[is.na(t21b)] = 0
t21c = cbind(t21b, pct_cov=t21b$qLen/t21b$len_scaf)
t21d = t21c[t21c$qLen<1000,]

f21 = file.path(dir, "21_scaffold_status.tbl")
write.table(t21c, file=f21, col.names=T, row.names=F, sep="\t", quote=F)
write.table(t21d, file=file.path(dir, "22_unmapped.tbl"), col.names=T, row.names=F, sep="\t", quote=F)

##### run assembly.pl to blast NR/NT database
f36 = file.path(dir, "36_annotated.tbl")
t36 = read.table(f36, header=TRUE, sep="\t", as.is=T)
t37 = ddply(t36, .(qId), summarise, qLen=sum(qLen), hLen=sum(hLen), n_seg=length(hId), pct=mean(pct), e=min(e), score = sum(score), n_tax=length(unique(category)), category=category[which(score==max(score))[1]], species=species[1])

t41 = cbind(t21c, category='', stringsAsFactors=FALSE)
t41$category[t21c$qLen>0] <- 'Medicago'
for (i in 1:nrow(t37)) {
  j = which(t41$qId == t37$qId[i])
  if(t41$score[j] < t37$score[i]) {
    t41[j, 2:7] = t37[i, 2:7]
    t41$category[j] = t37$category[i]
    t41$pct_cov[j] = t41$qLen[j] / t41$len_scaf[j]
  }
}
table(t41$category)

scaffold_stat <- function(ids, df) {
  dfs = df[df$qId %in% ids,]
  c('n'=length(ids), 'len_scaf'=sum(dfs$len_scaf), 'qLen_aln'=sum(dfs$qLen), 'hLen_aln'=sum(dfs$hLen))
}
id_all = t41$qId
id_mt = t41$qId[t41$category == "Medicago"]
id_fa = t41$qId[t41$category == "Fabaceae"]
id_un = t41$qId[t41$category == ""]
id_misc = id_all[! id_all %in% c(id_mt, id_fa, id_un)]
sapply(list(all=id_all, medicago=id_mt, fabaceae=id_fa, unknown=id_un, misc=id_misc), scaffold_stat, t41)

f41 = file.path(dir, "41_scaffold_status.tbl")
write.table(t41, file=f41, col.names=T, row.names=F, sep="\t", quote=F)
t41b = t41[t41$qId %in% c(id_fa, id_un),]
write.table(t41b, file=file.path(dir, "42_unknown.tbl"), col.names=T, row.names=F, sep="\t", quote=F)


