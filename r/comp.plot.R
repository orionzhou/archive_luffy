library(ggplot2)
library(plyr)
library(grid)
library(RColorBrewer)
source("comp.plot.fun.R")

f_gap_ref = file.path(DIR_Data, "genome/Mtruncatula_4.0/52_gap_loc.tbl")
t_gap_ref = read.table(f_gap_ref, header=TRUE, sep="\t", as.is=T)
#tg = tg[tg$len>=5000 & substr(tg$id,1,3)=="chr",]

f_len_ref = file.path(DIR_Data, "genome/Mtruncatula_4.0/51_seqlen.tbl")
t_len_ref = read.table(f_len_ref, header=TRUE, sep="\t", as.is=T)

f_crp = file.path(DIR_Misc4, "spada.crp.Mtruncatula_4.0/31_model_evaluation/61_final.gtb")
t_crp = read.table(f_crp, header=TRUE, sep="\t", as.is=T)[,c(1:6,13:18)]

f_gen_ref = file.path(DIR_Data, "genome/Mtruncatula_4.0/65_phase_fixed.gtb")
t_gen_ref = read.table(f_gen_ref, header=TRUE, sep="\t", as.is=T, quote=NULL)[,c(1:6,13:18)]


acc = "hm056"
#acc = "hm340"
dir = file.path(DIR_Misc3, acc)

f_len = file.path(dir, "11_seqlen.tbl")
t_len = read.table(f_len, header=TRUE, sep="\t", as.is=T)

f_gap = file.path(dir, "12_gaploc.tbl")
t_gap = read.table(f_gap, header=TRUE, sep="\t", as.is=T)


f05 = file.path(dir, '21_blastn/05_tiled.tbl')
t05 = read.table(f05, header=TRUE, sep="\t", as.is=T)
t_aln = t05

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

##### plot a chromosome region with mapping info
title = sprintf("%s / Mt4.0 (A17)", toupper(acc))

source("assembly_functions.R")
chr = 'chr7'
beg = 20022555
end = 20116342
tas = t_aln[t_aln$hId==chr & ( (beg<=t_aln$hBeg & t_aln$hBeg<=end) | (beg<=t_aln$hEnd & t_aln$hEnd<=end) ), ]
dat = data_preprocess(tas, t_len, t_len_ref, t_gap, t_gap_ref, t_gen_ref)
fn = sprintf("%s/figs/%s_%d_%d.png", dir, chr, beg, end)
plot_final(fn, dat, width=2000, height=1000, title=title)

source("assembly_functions.R")
nums = c(0,1,119,189,359,432,1000,1180,1181,1186,1188,1191,1203,1824,3189,4144)
#nums = c(1203)
for (num in nums) {
	id = sprintf("scaffold_%d", num)
	fn = sprintf("%s/figs/%s.png", dir, id)
	tas = t_aln[t_aln$qId==id & t_aln$qLen >= 1000,]
	dat = data_preprocess(tas, t_len, t_len_ref, t_gap, t_gap_ref, t_gen_ref)
	plot_final(fn, dat, width=2000, height=1000, title=title)
}
