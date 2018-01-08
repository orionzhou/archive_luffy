require(plyr)
require(ggplot2)
require(dplyr)
require(tidyr)
require(Hmisc)
options(scipen=999)

# process B73 TE GFF
dirg = '/home/springer/zhoux379/data/genome/B73/TE'
fi = file.path(dirg, "01.gff")
ti = read.table(fi, sep = "\t", as.is = T, header = F, quote = '')
colnames(ti) = c("chr", "src", "type", "beg", "end", "score", "srd", "phase", "note")
tmp = sapply(strsplit(ti$note, split = ";", fixed = T), "[", 1)
tids = sapply(strsplit(tmp, split = "=", fixed = T), "[", 2)
stopifnot(length(tids)==length(unique(tids)))

to = cbind(ti[,c(1,4,5)], tid = tids)
fo = sprintf("%s/05.tsv", dirg)
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

# process W22 TE GFF
dirg = '/home/springer/zhoux379/data/genome/W22/TE'
fi = file.path(dirg, "01.gff")
ti = read.table(fi, sep = "\t", as.is = T, header = F, quote = '')
colnames(ti) = c("chr", "src", "type", "beg", "end", "score", "srd", "phase", "note")
tids = sapply(strsplit(ti$note, split = "=", fixed = T), "[", 2)
stopifnot(length(tids)==length(unique(tids)))

ti$chr[ti$chr == 'chr10000001'] = 'unmapped'
to = cbind(ti[,c(1,4,5)], tid = tids)
fo = sprintf("%s/05.tsv", dirg)
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)


# generate flanking sequences
dirw = '/home/springer/zhoux379/data/misc1/te.flanking'
ft = file.path(dirw, "00.cfg.tsv")
tt = read.table(ft, sep = "\t", as.is = T, header = T)

tt1 = unique(tt[,-2])
for (i in 1:nrow(tt1)) {
    opt = tt1$opt[i]
    gt = tt1$genotype[i]
    tgt = tt1$tgt[i]
    
    fl = sprintf('/home/springer/zhoux379/data/genome/%s/TE/05.tsv', gt)
	tl = read.table(fl, sep = "\t", as.is = T, header = T)

    tb1 = tl
    tb1$beg = tl$beg + tt$u_beg[i]
    tb1$end = tl$beg + tt$u_end[i]
    tb1$beg = tb1$beg - 1
    
    tb2 = tl
    tb2$beg = tl$end + tt$d_beg[i]
    tb2$end = tl$end + tt$d_end[i]
    tb2$beg = tb2$beg - 1
    
    tb = rbind(tb1, tb2)
    tb = tb[order(tb$tid),]
    sum(tb$beg < 0)
    
    fb1 = sprintf("%s/01.loc/%s.%s.1.bed", dirw, gt, opt)
    write.table(tb1, fb1, sep = "\t", row.names = F, col.names = F, quote = F)
    fb2 = sprintf("%s/01.loc/%s.%s.2.bed", dirw, gt, opt)
    write.table(tb2, fb2, sep = "\t", row.names = F, col.names = F, quote = F)
    
    pre = sprintf("%s.%s", gt, opt)
    cmd = sprintf("seqret.py --padding %s %s %s/02.seq/%s.1.fas", gt, fb1, dirw, pre)
    cat(cmd, "\n")
    system(cmd)
    cmd = sprintf("seqret.py --padding %s %s %s/02.seq/%s.2.fas", gt, fb2, dirw, pre)
    cat(cmd, "\n")
    system(cmd)
}


for (i in 1:nrow(tt)) {
    opt = tt$opt[i]
    gt = tt$genotype[i]
    tgt = tt$tgt[i]
    
    pre1 = sprintf("%s.%s", gt, opt)
    pre2 = sprintf("%s.%s.%s", gt, tgt, opt)
    cmd = sprintf("bwa mem -t 8 $genome/%s/21.bwa/db %s/02.seq/%s.1.fas %s/02.seq/%s.2.fas > %s/04.bwa/%s.sam", tgt, dirw, pre1, dirw, pre1, dirw, pre2)
    #cat(cmd, "\n")
    
    cmd = sprintf("sam2tsv.py %s/04.bwa/%s.sam %s/04.bwa/%s.tsv", dirw, pre2, dirw, pre2)
    cat(cmd, "\n")
    system(cmd)
    cmd = sprintf("psl.filter.py --ident 0.9 --cov 0.9 --best %s/04.bwa/%s.tsv %s/07.filtered/%s.tsv", dirw, pre2, dirw, pre2)
    cat(cmd, "\n")
    system(cmd)
}
# cd $misc1/te.flanking
# seq.splitn.py --n 8 --cpu 8 02.fas 04.blat.W22
# seq.splitn.py --n 8 --cpu 8 02.fas 04.blat.PH207
# cat 04.blat.W22/*.psl > 05.blat.W22.psl
# cat 04.blat.PH207/*.psl > 05.blat.PH207.psl
# psl2tsv.pl -i 05.blat.W22.psl -o 06.W22.tsv
# psl2tsv.pl -i 05.blat.PH207.psl -o 06.PH207.tsv
# psl.filter.py --ident 0.9 --cov 0.9 --best 06.W22.tsv 07.W22.tsv
# psl.filter.py --ident 0.9 --cov 0.9 --best 06.PH207.tsv 07.PH207.tsv


i = 1
    opt = tt$opt[i]
    gt = tt$genotype[i]
    tgt = tt$tgt[i]

    fl = sprintf('/home/springer/zhoux379/data/genome/%s/TE/05.tsv', gt)
	tl = read.table(fl, sep = "\t", as.is = T, header = T)
	tl = cbind(tl, rsize = tl$end - tl$beg + 1)
	ids_all = tl$tid

fi = sprintf("%s/07.filtered/%s.%s.%s.tsv", dirw, gt, tgt, opt)
tir = read.table(fi, sep = "\t", header = T, as.is = T)[,1:20]
ti = tir[tir$alnLen/tir$qSize >= 0.9 & tir$ident >= 0.9,]

res = strsplit(ti$qId, split = ".", fixed = T)
tids = sapply(res, "[", 1)
pe = sapply(res, "[", 2)
ti = cbind(ti, tid = tids, pe = sprintf("flank%s", pe))

ids_mapped = unique(ti$qId)

grp = dplyr::group_by(ti, qId)
ti2 = dplyr::summarise(grp, nbest = sum(score == max(score)), score = max(score))
stopifnot(sum(ti2$nbest) == nrow(ti))
ti3 = merge(ti, ti2, by = c("qId", "score"))
ids_mapped_u = unique(ti2$qId[ti2$nbest == 1])
ids_mapped_m = unique(ti2$qId[ti2$nbest > 1])

# unique hits
tk = ti3[ti3$qId %in% ids_mapped_u,]
stopifnot(nrow(tk) == length(ids_mapped_u))
ids_mapped_u1 = tk$qId[tk$alnLen == tk$qSize & tk$misMatch == 0 & tk$qNumIns+tk$tNumIns==0]
ids_mapped_u2 = tk$qId[tk$alnLen < tk$qSize | tk$misMatch > 0 | tk$qNumIns+tk$tNumIns > 0]
stopifnot(length(ids_mapped_u1) + length(ids_mapped_u2) == length(ids_mapped_u))


outputs = c(
'',
sprintf("%6d total %s TE sequences:", length(ids_all), gt),
sprintf("%6d do not map to %s (0.9 identity, 0.9 coverage)", tgt, length(ids_all)-length(ids_mapped)),
sprintf("%6d mapped >=1 times:", length(ids_mapped)),
sprintf("   %6d mapped uniquely:", length(ids_mapped_u)),
sprintf("      %6d are identical:", length(ids_mapped_u1)),
sprintf("      %6d have at least 1 mismatch/indel:", length(ids_mapped_u2)),
sprintf("   %6d mapped multiple times:", length(ids_mapped_m)),
''
)
cat(paste(outputs, collapse = "\n"))


## look at pair mapping
grp = dplyr::group_by(ti, tid, pe)
tb = dplyr::summarise(grp, nbest = sum(score == max(score)), score = max(score))
stopifnot(sum(tb$nbest) == nrow(ti))

tp = spread(tb[,c(1:4)], pe, nbest)
tp$flank1[is.na(tp$flank1)] = 0
tp$flank2[is.na(tp$flank2)] = 0

ids_mapped = unique(tp$tid)
n_unmap = length(ids_all) - length(ids_mapped)

tp1 = tp[tp$flank1 == 0 | tp$flank2 == 0,]
tp2 = tp[tp$flank1 > 0 & tp$flank2 > 0,]
stopifnot(nrow(tp1)+nrow(tp2) == nrow(tp))

tp11 = tp1[tp1$flank1 + tp1$flank2 == 1,]
tp12 = tp1[tp1$flank1 + tp1$flank2 > 1,]
stopifnot(nrow(tp11)+nrow(tp12) == nrow(tp1))

tp21 = tp2[tp2$flank1 == 1 & tp2$flank2 == 1,]
tp22 = tp2[(tp2$flank1 == 1 & tp2$flank2 > 1) | (tp2$flank1 > 1 & tp2$flank2 == 1),]
tp23 = tp2[tp2$flank1 > 1 & tp2$flank2 > 1,]
stopifnot(nrow(tp21)+nrow(tp22)+nrow(tp23) == nrow(tp2))


# look at mapped flank distance
td = ti[ti$tid %in% tp21$tid,]
stopifnot(nrow(td) == 2*nrow(tp21))
#td$rsize = as.numeric(as.character(td$rsize))
td1 = td[td$pe == 'flank1',]
td2 = td[td$pe == 'flank2',]
td1 = td1[order(td1$tid),]
td2 = td2[order(td2$tid),]
identical(td1$tid, td2$tid)
td3 = cbind(td1[,-1], schr = td1$tId==td2$tId, asize = pmin(abs(td1$tBeg-td2$tEnd-1), abs(td1$tEnd-td2$tBeg-1)))

td4 = merge(td3, tl[,c("tid","rsize")], by = 'tid')
td4 = td4[td4$schr,]
describe(td4$rsize-td4$asize)
dists = abs(td4$rsize-td4$asize)
sum(dists > 10000)
td5 = td4[td4$asize <= 100000,]

p1 = ggplot(td5) +
  geom_point(aes(x = rsize, y = asize), size = 1) +
  scale_x_continuous(name = 'Gap size in B73 genome') +
  scale_y_continuous(name = sprintf("Gap size in %s genome", gt)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 10), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8))
fp = sprintf("%s/11.%s.pdf", dirw, gt)
ggsave(p1, filename = fp, width = 8, height = 8)


fo = sprintf("%s/10.%s.tsv", dirw, gt)
write.table(td5[,21:24], fo, sep = "\t", row.names = F, col.names = T, quote = F)

sum(td5$asize < 10)



outputs = c(
'',
sprintf("%6d total TEs:", length(ids_all)),
sprintf("%6d: neither flank mapped to %s (0.9 identity, 0.9 coverage)", n_unmap, gt),
sprintf("%6d: only one flank mapped", nrow(tp1)),
sprintf("   %6d mapped uniquely", nrow(tp11)),
sprintf("   %6d mapped multiple times", nrow(tp12)),
sprintf("%6d: both flanks mapped", nrow(tp2)),
sprintf("   %6d both flanks mapped uniquely", nrow(tp21)),
sprintf("     %6d uniquely mapped to the same chromosome", nrow(td4)),
sprintf("       %d with gap size < 100kb => used to make the attached figure", nrow(td5)),
sprintf("   %6d one flank mapped uniquely but the other mapped multiple times", nrow(tp22)),
sprintf("   %6d both flanks mapped multiple times", nrow(tp23)),
''
)
cat(paste(outputs, collapse = "\n"))
