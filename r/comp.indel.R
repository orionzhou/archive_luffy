require(rtracklayer)
source("comp.fun.R")

tname = "hm101"
qname1 = "hm056"
qname2 = "hm340"
qname = qname1

t = read_genome_stat(tname)
q = read_genome_stat(qname)
cq = read_comp_stat(qname, tname)

dir = sprintf("%s/%s", cq$dir, "23_blat")

# filter SVs overlapping with qGap/tGap
tidm = read.table(file.path(dir, "26.vnt/idm"), header = F, sep = "\t", as.is = T)
colnames(tidm) = c('tId', 'tBeg', 'tEnd', 'id', 'qId', 'qBeg', 'qEnd')
tidm = cbind(tidm, qlen = tidm$qEnd - tidm$qBeg - 1, 
  tlen = tidm$tEnd - tidm$tBeg - 1)

rm_ovlp_t <- function(ds, gr) {
  grl = with(ds, makeGRangesListFromFeatureFragments(
    seqnames = tId, fragmentStarts = sprintf("%d,", tBeg + 1), 
    fragmentWidths = sprintf("%d,", tEnd - tBeg), 
    strand = rep("+", nrow(ds))))
  ma = as.matrix(findOverlaps(gr, grl))
  didx1 = data.frame(sidx=ma[,2], qidx=ma[,1])
  idxs_rm = unique(didx1$sidx)
  ds[-idxs_rm, ]
}
rm_ovlp_q <- function(ds, gr) {
  grl = with(ds, makeGRangesListFromFeatureFragments(
    seqnames = qId, fragmentStarts = sprintf("%d,", qBeg + 1), 
    fragmentWidths = sprintf("%d,", qEnd - qBeg), 
    strand = rep("+", nrow(ds))))
  ma = as.matrix(findOverlaps(gr, grl))
  didx1 = data.frame(sidx=ma[,2], qidx=ma[,1])
  idxs_rm = unique(didx1$sidx)
  ds[-idxs_rm, ]
}

gr_tgap = GRanges(seqnames = t$gap$id, 
  ranges = IRanges(t$gap$beg, end = t$gap$end))
gr_qgap = GRanges(seqnames = q$gap$id, 
  ranges = IRanges(q$gap$beg, end = q$gap$end))

tm = rm_ovlp_q(rm_ovlp_t(tidm, gr_tgap), gr_qgap)
#tm = tm[order(tm$tid, tm$tbeg, tm$tend), ]

nrow(tidm)
nrow(tm)

# filter small indels
tms = tm[tm$tlen < 100 & tm$qlen < 100, ]
tml = tm[tm$tlen >= 100 | tm$qlen >= 100, ]
nrow(tms)
nrow(tml)

write.table(tms, file.path(dir, "27.vnt/03.s.tbl"), sep = "\t",
  quote = F, row.names = F, col.names = T)
write.table(tml, file.path(dir, "27.vnt/03.l.tbl"), sep = "\t",
  quote = F, row.names = F, col.names = T)

# filter flanking regions & complex rearrangements
tm = read.table(file.path(dir, "27.vnt/03.l.tbl"), header = T, sep = "\t")

i = 1

tid = tm$tId[i]
tbeg = tm$tBeg[i]
tend = tm$tEnd[i]
qid = tm$qId[i]
qbeg = tm$qBeg[i]
qend = tm$qEnd[i]



if(tend - tbeg - 1 > 0) {
  tt = read_gax(cq$tgax, cq$tsnp, tid, tbeg, tend, 't')
  qlb = qbeg - 500 + 1
  qle = qbeg
  qrb = qend
  qre = qend + 500 - 1
}
if(qend - qbeg - 1 > 0) {
  tq = read_gax(cq$qgax, cq$qsnp, qid, qbeg, qend, 'q')
}



# plot
ts = cbind(tc, lenlog=log(tc$aLen))
p <- ggplot(ts) +
  geom_point(mapping=aes(x=qCov, y=hCov, color=lenlog), size=1) +
  scale_colour_gradient(low="white", high="red")

qid.levels = unique(tc$qId[order(tc$hId, tc$hBeg)])
dat1.coord = data.frame(id=qid.levels, order=1:length(qid.levels))
dat1.coord = merge(dat1.coord, dat1.seqlen, by='id')
dat1.coord = dat1.coord[order(dat1.coord$order),-2]
dat1.coord = cbind(dat1.coord, beg=1, end=dat1.coord$length)
for (i in 2:nrow(dat1.coord)) {
  dat1.coord$beg[i] = dat1.coord$end[i-1] + 1
  dat1.coord$end[i] = dat1.coord$beg[i] + dat1.coord$length[i] - 1
}

colnames(dat1.coord) = c('qId', 'qLen.r', 'qBeg.r', 'qEnd.r')
alnplot = merge(alnplot, dat1.coord, by='qId')
colnames(dat2.coord) = c('hId', 'hLen.r', 'hBeg.r', 'hEnd.r')
alnplot = merge(alnplot, dat2.coord, by='hId')

p <- ggplot(tc) +
  geom_segment(mapping=aes(x=hBeg, xend=hEnd, y=qBeg, yend=qEnd), size=0.6) +
  layer(data=dat2.coord, geom='rect', mapping=aes(xmin=hBeg.r, xmax=hEnd.r, 
    ymin=-3000000, ymax=0, fill=hId), geom_params=list(size=0)) +
  layer(data=dat2.coord, geom='text', mapping=aes(
    x=(hBeg.r+hEnd.r)/2, y=-4000000, label=hId), 
    geom_params=list(hjust=0.5, vjust=1, angle=0, size=3)) +
  facet_grid(qId ~ hId, scales="free") +
  scale_color_brewer(palette="Set1") +
  scale_x_continuous(name='HM101 (Mt4.0)') +
  scale_y_continuous(name=dat1.name) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank())
ggsave(file.path(dat1.dir, "all.png"), p, width=8, height=8)



f23_21 = file.path(dat1.dir, "23_blat/21_indel.gal")
tid1 = read.table(f23_21, header=TRUE, sep="\t", as.is=T)

tid = tid1
tail(sort(tid$hLen), 100)


tid2 = ddply(tid1, .(idc), summarise, qId=unique(qId), qBeg=min(qBeg), qEnd=max(qEnd), qSrd=unique(qSrd), qLen=sum(qLen), hId=unique(hId), hBeg=min(hBeg), hEnd=max(hEnd), hSrd=unique(hSrd), hLen=sum(hLen))
tid3 = cbind(tid2, qCov=tid2$qLen/(tid2$qEnd-tid2$qBeg+1), hCov=tid2$hLen/(tid2$hEnd-tid2$hBeg+1))

##### identify insertion / deletions
indel = data.frame()
gap1 = t_gap
gap2 = t_gap_ref
for (i in 2:nrow(t_aln)) {
  df = t_aln
  if( df$qId[i]==df$qId[i-1] & df$hId[i]==df$hId[i-1] & df$strand[i]==df$strand[i-1] ) {
    qId = df$qId[i]; hId = df$hId[i]; strand = df$strand[i]
    qBeg1 = df$qBeg[i-1]; qEnd1 = df$qEnd[i-1]; qBeg2 = df$qBeg[i]; qEnd2 = df$qEnd[i]
    hBeg1 = df$hBeg[i-1]; hEnd1 = df$hEnd[i-1]; hBeg2 = df$hBeg[i]; hEnd2 = df$hEnd[i]
    if( strand == "-" ) {
      qItv = qBeg2 - qEnd1 - 1
      hItv = hBeg1 - hEnd2 - 1
      qHasGap = sum( gap1$id==qId & gap1$beg >= qEnd1 & gap1$end <= qBeg2 )
      hHasGap = sum( gap2$id==hId & gap2$beg >= hEnd2 & gap2$end <= hBeg1 )
    } else {
      qItv = qBeg2 - qEnd1 - 1
      hItv = hBeg2 - hEnd1 - 1
      qHasGap = sum( gap1$id==qId & gap1$beg >= qEnd1 & gap1$end <= qBeg2 )
      hHasGap = sum( gap2$id==hId & gap2$beg >= hEnd1 & gap2$end <= hBeg2 )
    }
    type = ''
    if( abs(qItv) <= 10 & hItv >= 100 & hItv <= 10000 & hHasGap == 0) {
      type = 'del'
    } else if( abs(hItv) <= 10 & qItv >= 100 & qItv <= 10000 & qHasGap == 0) {
      type = "ins"
    }
    if(type != '') {
      indel.1 = data.frame('qId'=qId, 'qBeg1'=qBeg1, 'qEnd1'=qEnd1, 'qBeg2'=qBeg2, 'qEnd2'=qEnd2, 'strand'=strand, 
      'hId'=hId, 'hBeg1'=hBeg1, 'hEnd1'=hEnd1, 'hBeg2'=hBeg2, 'hEnd2'=hEnd2, 
      'qItv'=qItv, 'hItv'=hItv, 'qLen1'=df$qLen[i-1], 'qLen2'=df$qLen[i], 
      'hLen1'=df$hLen[i-1], 'hLen2'=df$hLen[i], 'pct1'=df$pct[i-1], 'pct2'=df$pct[i], 
      'score1'=df$score[i-1], 'score2'=df$score[i], 'type'=type)
      indel = rbind(indel, indel.1)
    }
  }
}
nrow(indel)
head(indel)

f31_01 = file.path(dir, "31_indel/01_raw.tbl")
write.table(indel, file=f31_01, col.names=T, row.names=F, sep="\t", quote=F)

##### read indels
dat1.name = "hm056"
dat1.dir = file.path("/home/youngn/zhoup/Data/misc3", dat1.name)

f31_01 = file.path(dat1.dir, "31_indel/01_raw.tbl")
indel = read.table(f31_01, header=TRUE, sep="\t", as.is=T)

idxs = c()
dfq = indel
tg = rbind(dat2$crp, dat2$nbs)
for( i in 1:nrow(dfq) ) {
  chr = dfq$hId[i];  beg = dfq$hEnd1[i]; end = dfq$hBeg2[i]
  flags = tg$chr==chr & ( (beg < tg$beg & tg$beg < end) | (beg < tg$end & tg$end < end) )
  if(sum(flags) > 0) {
    idxs = c(idxs, i)
  }
}
indel.1 = indel[idxs,]
nrow(indel.1)
head(indel.1)

dfp = indel.1
for (i in 1:nrow(dfp)) {
  chr = dfp$hId[i]
  beg = min( dfp$hBeg1[i], dfp$hBeg2[i] )
  end = max( dfp$hEnd1[i], dfp$hEnd2[i] )
  tas = aln[aln$hId==chr & ((beg<=aln$hBeg & aln$hBeg<=end) | (beg<=aln$hEnd & aln$hEnd<=end)), ]
  dat = data_preprocess(tas, dat1, dat2)
  fn = sprintf("%s/figs/indel_%d.png", dat1.dir, i)
  subtitle = sprintf("%s:%g-%gK", chr, beg/1000, end/1000)
  comp.plot(fn, dat, width=2000, height=1000, subtitle=subtitle)
}


idxs = c()
dir_pindel = "/home/youngn/zhoup/Data/misc3/hapmap_mt40/40_sv"
f_pindel_del = file.path(dir_pindel, "HM056_chr5_D.tbl")
td = read.table(f_pindel_del, head=T, sep="\t", as.is=T)
td.1 = td
td.2 = td.1[abs(td.1$svlen) >= 1 & abs(td.1$svlen) <= 10000, 1:7]

idxs = c()
dfq = td.2
tg = dat2$gene
tg = rbind(dat2$crp, dat2$nbs)
for( i in 1:nrow(dfq) ) {
  chr = dfq$chr[i]; beg = dfq$pos[i]; end =  dfq$end[i]
  flags = tg$chr==chr & ( (beg < tg$beg & tg$beg < end) | (beg < tg$end & tg$end < end) )
  if(sum(flags) > 0) { idxs = c(idxs, i) }
}
td.3 = dfq[idxs,]
nrow(td.3)
head(td.3)
