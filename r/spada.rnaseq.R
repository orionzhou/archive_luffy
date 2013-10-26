library(ggplot2)
library(plyr)
library(data.table)

org = "Athaliana"
org = "Mtruncatula_3.5"

dir = file.path(DIR_Misc2, "spada.rnaseq", org, "31_rnaseq")

fi1 = file.path(dir, "05_gs.tbl")
fi2 = file.path(DIR_Misc2, "spada.affy", org, "25.tbl")
fo = file.path(dir, "06_gs.tbl")

fi1 = file.path(dir, "15_spada.tbl")
fi2 = file.path(DIR_Misc2, "spada.affy", org, "35.tbl")
fo = file.path(dir, "16_spada.tbl")

tr = read.table(fi1, sep="\t", header=T, as.is=T)[,-4]
tr$fpkm[is.na(tr$fpkm)] = 0
tr$cov[is.na(tr$cov)] = 0
#tr = cbind(tr, call.rna=tr$fpkm>=1)
tr = cbind(tr, call.rna=tr$cov>=2)

ta = read.table(fi2, sep="\t", header=T, as.is=T)
identical(tr$id, ta$id)
t = cbind(tr, ta[-1])
t = cbind(t, call=t$call.rna | t$call.affy)
sum(t$call, na.rm=T)
sum(t$call, na.rm=T)/nrow(t)

write.table(t, fo, row.names=F, col.names=T, sep="\t", quote=F)
