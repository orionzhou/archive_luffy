library("GenomicRanges")
library("plyr")

# load query stats
org.q = "HM056"
dir.q = file.path('/home/youngn/zhoup/Data/genome', org.q)

fs.q = file.path(dir.q, '15_seqlen.tbl')
fg.q = file.path(dir.q, '16_gaploc.tbl')

ts.q = read.table(fs.q, sep='\t', header=T, as.is=T)
tg.q = read.table(fg.q, sep='\t', header=T, as.is=T)

si.q = Seqinfo(ts.q$id, seqlengths=ts.q$length)
gra.q = GRanges(seqnames=Rle(ts.q$id), ranges=IRanges(rep(1, nrow(ts.q)), end=ts.q$length), seqinfo=si.q)
grg.q = GRanges(seqnames=Rle(tg.q$id), ranges=IRanges(tg.q$beg, end=tg.q$end), seqinfo=si.q)

sum(width(gra.q))
sum(width(grg.q))

# load target stats
org.t = "HM101"
dir.t = file.path('/home/youngn/zhoup/Data/genome', org.t)

fs.t = file.path(dir.t, '15_seqlen.tbl')
fg.t = file.path(dir.t, '16_gaploc.tbl')

ts.t = read.table(fs.t, sep='\t', header=T, as.is=T)
tg.t = read.table(fg.t, sep='\t', header=T, as.is=T)

si.t = Seqinfo(ts.t$id, seqlengths=ts.t$length)
gra.t = GRanges(seqnames=Rle(ts.t$id), ranges=IRanges(rep(1, nrow(ts.t)), end=ts.t$length), seqinfo=si.t)
grg.t = GRanges(seqnames=Rle(tg.t$id), ranges=IRanges(tg.t$beg, end=tg.t$end), seqinfo=si.t)

sum(width(gra.t))
sum(width(grg.t))

# load comparison
dir = sprintf('/home/youngn/zhoup/Data/misc3/%s_%s', org.q, org.t)
# round 1 - generating novel sequences
# galexpand.pl -i 35.gal -o 35.gall
t1w = read.table(file.path(dir, '23_blat/35.gal'), sep='\t', header=T, as.is=T)[,1:17]
sum(t1w$match) / (sum(t1w$match) + sum(t1w$misMatch))
t1 = read.table(file.path(dir, '23_blat/35.gall'), sep='\t', header=T, as.is=T)

gr1.1 = GRanges(seqnames=Rle(t1$tId), ranges=IRanges(t1$tBeg, end=t1$tEnd), seqinfo=si.q)
gr1.2 = reduce(gr1.1)
sum(width(gr1.2))

grn1 = setdiff(gra.q, union(gr1.2, grg.q))
sum(width(grn1))

tn1.1 = data.frame(id=as.character(seqnames(grn1)), beg=as.numeric(start(grn1)), end=as.numeric(end(grn1)), len=as.numeric(width(grn1)))
tn1.2 = tn1.1[tn1.1$len>=1000,]
sum(tn1.2$len)
write.table(tn1.2, file.path(dir, '41_novelseq/10_novel.tbl'), col.names=T, row.names=F, sep='\t', quote=F)


# round 2
# seqextract.pl -i $data/genome/$org_q/11_genome.fa -o 11_novel.fa -n 10_novel.tbl
# blat $data/db/blat/Mtruncatula_4.0.2bit -ooc=$data/db/blat/Mtruncatula_4.0.2bit.tile11.ooc 11_novel.fa 12.psl -noHead -noTrimA
# psl2gal.pl -i 12.psl -o - | galcoord.pl -i - -p qry -o - | galfix.pl -i - -q $data/genome/$org_q/11_genome.fa -t $data/genome/HM101/11_genome.fa -o 13.gal
# galexpand.pl -i 13.gal -o 13.gall
t2w = read.table(file.path(dir, '41_novelseq/13.gal'), sep='\t', header=T, as.is=T)[,1:17]
sum(t2w$match) / (sum(t2w$match) + sum(t2w$misMatch))
t2 = read.table(file.path(dir, '41_novelseq/13.gall'), sep='\t', header=T, as.is=T)

gr2.1 = GRanges(seqnames=Rle(t2$qId), ranges=IRanges(t2$qBeg, end=t2$qEnd), seqinfo=si.q)
gr2.2 = reduce(gr2.1)
sum(width(gr2.2))

grn2 = setdiff(grn1, gr2.2)
sum(width(grn2))

tn2.1 = data.frame(id=as.character(seqnames(grn2)), beg=as.numeric(start(grn2)), end=as.numeric(end(grn2)), len=as.numeric(width(grn2)))
tn2.2 = tn2.1[tn2.1$len>=1000,]
sum(tn2.2$len)
# write.table(tn2.2, file.path(dir, '41_novelseq/20_novel.tbl'), col.names=T, row.names=F, sep='\t', quote=F)

# round 3
# seqextract.pl -i $data/genome/$org_q/11_genome.fa -o 21_novel.fa -n 20_novel.tbl
# blat $data/db/blat/Mtruncatula_4.0.2bit 21_novel.fa -ooc=$data/db/blat/Mtruncatula_4.0.2bit.tile11.ooc 22.psl -noHead -noTrimA
# psl2gal.pl -i 22.psl -o - | galcoord.pl -i - -p qry -o - | galfix.pl -i - -q $data/genome/$org_q/11_genome.fa -t $data/genome/HM101/11_genome.fa -o 23.gal
# galexpand.pl -i 23.gal -o 23.gall
t3w = read.table(file.path(dir, '41_novelseq/23.gal'), sep='\t', header=T, as.is=T)[,1:17]
sum(t3w$match) / (sum(t3w$match) + sum(t3w$misMatch))
t3 = read.table(file.path(dir, '41_novelseq/23.gall'), sep='\t', header=T, as.is=T)

gr3.1 = GRanges(seqnames=Rle(t3$qId), ranges=IRanges(t3$qBeg, end=t3$qEnd), seqinfo=si.q)
gr3.2 = reduce(gr3.1)
sum(width(gr3.2))

grn3 = setdiff(grn2, gr3.2)
sum(width(grn3))

tn3.1 = data.frame(id=as.character(seqnames(grn3)), beg=as.numeric(start(grn3)), end=as.numeric(end(grn3)), len=as.numeric(width(grn3)))
tn3.2 = tn3.1[tn3.1$len>=1000,]
sum(tn3.2$len)
# write.table(tn3.2, file.path(dir, '41_novelseq/30_novel.tbl'), col.names=T, row.names=F, sep='\t', quote=F)

# seqextract.pl -i $data/genome/$org_q/11_genome.fa -o 31_novel.fa -n 30_novel.tbl



# blastn validation
# blastn -db $data/db/blast/Mtruncatula_4.0 -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen length nident mismatch gaps evalue bitscore qseq sseq' -num_threads 4 -query 11_novel.fa -out 71_blastn.tbl
# blast2gal.pl -i 71_blastn.tbl -o - | galcoord.pl -i - -p qry -o - | galexpand.pl -i - -o 72.gall
# gal.pl -i 72.gal -o 73.gal -opt coordq
tv1 = read.table(file.path(dir, '41_novelseq/73.gal'), sep='\t', header=T, as.is=T)
sum(tv1$match)/sum(tv1$match/(tv1$ident/100))
grv1.1 = GRanges(seqnames=Rle(tv1$qId), ranges=IRanges(tv1$qBeg, end=tv1$qEnd), seqinfo=si.q)
grv1.2 = reduce(grv1.1)
sum(width(grv1.2))

# blastn -db $data/db/blast/Mtruncatula_4.0 -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen length nident mismatch gaps evalue bitscore qseq sseq' -num_threads 4 -query 21_novel.fa -out 81_blastn.tbl
# blastToGal.pl -i 81_blastn.tbl -o 82.gal
# gal.pl -i 82.gal -o 83.gal -opt coordq
tv2 = read.table(file.path(dir, '41_novelseq/83.gal'), sep='\t', header=T, as.is=T)
sum(tv2$match)/sum(tv2$match/(tv2$ident/100))
grv2.1 = GRanges(seqnames=Rle(tv2$qId), ranges=IRanges(tv2$qBeg, end=tv2$qEnd), seqinfo=si.q)
grv2.2 = reduce(grv2.1)
sum(width(grv2.2))

# blastn -db $data/db/blast/Mtruncatula_4.0 -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen length nident mismatch gaps evalue bitscore qseq sseq' -num_threads 4 -query 31_novel.fa -out 91_blastn.tbl
# blastToGal.pl -i 91_blastn.tbl -o 92.gal
# gal.pl -i 92.gal -o 93.gal -opt coordq
tv3 = read.table(file.path(dir, '41_novelseq/93.gal'), sep='\t', header=T, as.is=T)
sum(tv3$match)/sum(tv3$match/(tv3$ident/100))
grv3.1 = GRanges(seqnames=Rle(tv3$qId), ranges=IRanges(tv3$qBeg, end=tv3$qEnd), seqinfo=si.q)
grv3.2 = reduce(grv3.1)
sum(width(grv3.2))
tv3 = read.table(file.path(dir, '41_novelseq/94.gal'), sep='\t', header=T, as.is=T)



# blastn NT
dir = file.path("/home/youngn/zhoup/Data/misc3/pan3", org.q)
tg = read.table(file.path(dir, "15.gal"), header=T, sep="\t")
tc = read.table(file.path(dir, "16_cat.tbl"), header=T, sep="\t", as.is=T)
tg.1 = merge(tg, tc, by.x='tId', by.y='id')

tg.2 = ddply(tg.1, .(kingdom), summarise, len=sum(tLen))