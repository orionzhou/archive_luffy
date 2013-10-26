library(plyr)
library(data.table)

dir = file.path(DIR_Misc2, "spada.affy", "Athaliana")
te1 = read.table(file.path(dir, "../atdefl_pac.tbl"), sep="\t", header=T, as.is=T)
totalcall = apply(te1[,seq(3,61,by=2)], 1, function(x) "P" %in% x)
te = data.frame(set=te1$Probe.set, call=totalcall)

#fi = file.path(dir, "23.tbl")
fi = file.path(dir, "33.tbl")
#fo = file.path(dir, "25.tbl")
fo = file.path(dir, "35.tbl")
tg = read.table(fi, sep="\t", header=T, as.is=T)[,c('id','set','num')]
tt = merge(tg, te, by='set', all.x=T)
tt$call[is.na(tt$call)] = ''
tt = tt[order(tt$id),c(2,1,3,4)]
df = tt
colnames(df) = c("id", "probe.set", "probe.number", "call.affy")
write.table(df, fo, row.names=F, col.names=T, sep="\t", quote=F)



dir = file.path(DIR_Misc2, "spada.affy", "Mtruncatula_3.5")
te1 = read.table(file.path(dir, "../mtdefl_pac.tbl"), sep="\t", header=T, as.is=T)
totalcall = apply(te1[,seq(6,206,by=8)], 1, function(x) "P" %in% x)
te = data.frame(set=te1$Probe.Set, call=totalcall)

fi = file.path(dir, "23.tbl")
#fi = file.path(dir, "33.tbl")
fo = file.path(dir, "25.tbl")
#fo = file.path(dir, "35.tbl")
tg = read.table(fi, sep="\t", header=T, as.is=T)[,c('id','set','num')]
tt = merge(tg, te, by='set', all.x=T)
tt$call[is.na(tt$call)] = ''
tt = tt[order(tt$id),c(2,1,3,4)]
df = tt
colnames(df) = c("id", "probe.set", "probe.number", "call.affy")
write.table(df, fo, row.names=F, col.names=T, sep="\t", quote=F)