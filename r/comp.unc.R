require(rtracklayer)
source("comp.fun.R")

tname = "hm101"
qname = "hm056"

t = read_genome_stat(tname)
q = read_genome_stat(qname)

dir = '/home/youngn/zhoup/Data/misc3/HM056_HM101/23_blat'
f23 = file.path(dir, "23.gax")
f33 = file.path(dir, "33.gax")
t23 = read.table(f23, header = F, sep = "\t", as.is = T)
t33 = read.table(f33, header = F, sep = "\t", as.is = T)
colnames(t23) = c('tid','tbeg','tend','tsrd','id','qid','qbeg','qend','qsrd')
colnames(t33) = c('tid','tbeg','tend','tsrd','id','qid','qbeg','qend','qsrd')


seqinfo = Seqinfo(q$len$id, seqlengths=q$len$size)
grt = GRanges(seqnames = q$len$id, 
  ranges = IRanges(rep(1, nrow(q$len)), end = q$len$size), seqinfo = seqinfo)
grg = GRanges(seqnames = q$gap$id, 
  ranges = IRanges(q$gap$beg, end = q$gap$end), seqinfo = seqinfo)   
 
sum(width(grt)) 
sum(width(grg)) 
  
gr23 = GRanges(seqnames = t23$qid, 
  ranges = IRanges(t23$qbeg, end = t23$qend), seqinfo = seqinfo)
gr33 = GRanges(seqnames = t33$tid, 
  ranges = IRanges(t33$tbeg, end = t33$tend), seqinfo = seqinfo)
sum(width(reduce(gr23)))
sum(width(reduce(gr33)))
sum(width(reduce(intersect(gr23, gr33))))
nov = setdiff(grt, union(grg, reduce(gr33)))
novs = nov[width(nov) >= 100]
sum(width(nov))
sum(width(novs))

gr33s = reduce(gr33)
t33s = data.frame(id = seqnames(gr33s), beg = start(gr33s)-1, end = end(gr33s))
write.table(t33s, file.path(dir, "tmp.bed"), col.names=F, row.names=F, sep="\t", quote=F)


d2 = '/home/youngn/zhoup/Data/misc3/HM056_HM101/41_novelseq'
fn1 = file.path(d2, 'nov1.tbl')
fn2 = file.path(d2, 'nov2.tbl')
fn3 = file.path(d2, 'nov3.tbl')
tn1 = read.table(fn1, header = T, sep = "\t", as.is = T)
tn2 = read.table(fn2, header = T, sep = "\t", as.is = T)
tn3 = read.table(fn3, header = T, sep = "\t", as.is = T)
grn1 = GRanges(seqnames = tn1$id, 
  ranges = IRanges(tn1$beg, end = tn1$end), seqinfo = seqinfo)
grn2 = GRanges(seqnames = tn2$id, 
  ranges = IRanges(tn2$beg, end = tn2$end), seqinfo = seqinfo)
grn3 = GRanges(seqnames = tn3$id, 
  ranges = IRanges(tn3$beg, end = tn3$end), seqinfo = seqinfo)

sum(width(grn1))
sum(width(grn2))
sum(width(grn3))

sum(width(intersect(gr33, grn1)))
sum(width(intersect(gr33, grn2)))
sum(width(intersect(gr33, grn3)))