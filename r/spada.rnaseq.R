library(ggplot2)
library(plyr)
library(data.table)

org = "Athaliana"
org = "Mtruncatula_3.5"

dir = file.path(DIR_Misc2, "spada.rnaseq", org, "31_rnaseq")
tg = read.table(file.path(dir, "05_gs.tbl"), sep="\t", header=T, as.is=T)
ts = read.table(file.path(dir, "15_spada.tbl"), sep="\t", header=T, as.is=T)


t = ts
cutoff_fpkm = 1 
sum(t$fpkm>cutoff_fpkm, na.rm=T)
sum(t$fpkm>cutoff_fpkm, na.rm=T)/nrow(t)

cutoff_cov = 2
sum(t$cov>cutoff_cov, na.rm=T)
sum(t$cov>cutoff_cov, na.rm=T)/nrow(t)


