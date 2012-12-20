library(seqinr)


dir = file.path(DIR_Misc2, "crp.evo")
f01 = file.path(dir, "01_seq.tbl")
ts = read.table(f01, header=T, sep="\t", as.is=T)

f09 = file.path(dir, "09.tbl")
t = read.table(f09, header=T, sep="\t", as.is=T)
t1 = t[,c('id','fam','id1','chr1','beg1','end1','str1','id2','chr2','beg2','end2','str2')]
chrs = unique(c(t1$chr1, t1$chr2))
t1$chr1 = as.numeric(factor(t1$chr1, levels=chrs))
t1$chr2 = as.numeric(factor(t1$chr2, levels=chrs))

id = 1
id1 = t1$id1[t1$id==id]
id2 = t1$id2[t1$id==id]
seq1 = s2c(ts$seq_ext[ts$id==id1])
seq2 = s2c(ts$seq_ext[ts$id==id2])
dotPlot(seq1,seq2,wsize=3)

t2 = t1[abs(t1$chr1-t1$chr2)==1,]
t3 = t1[t1$chr1==2|t1$chr2==1,]
t3 = t1[t1$chr1==2|t1$chr2==2,]
t3 = t1[t1$chr1==2|t1$chr2==3,]
t3 = t1[t1$chr1==2|t1$chr2==4,]
t3 = t1[t1$chr1==2|t1$chr2==5,]

ts = t3
xmax = max(ts$end1, ts$end2)
ymin = min(ts$chr1, ts$chr2)
ymax = max(ts$chr1, ts$chr2)
plot(0, type="n", xlim=c(0,xmax), ylim=c(ymin,ymax), xlab="", ylab="position", main="CRP duplication")
segments((ts$beg1+ts$end1)/2, as.numeric(ts$chr1), (ts$beg2+ts$end2)/2, as.numeric(ts$chr2))

