qorg = "HM056"
torg = "HM101"
dir = sprintf("%s/%s_%s/23_blat", DIR_Misc3, qorg, torg)

fg = sprintf("%s/genome/%s/51.gtb", DIR_Data, torg)
tg = read.table(fg, header=T, sep="\t", as.is=T, quote="")[,c(1,3:6,17)]
tg = cbind(tg, crp=grepl('CRP', tg$cat3), nbs=grepl('NBS', tg$cat3), te=grepl('TE', tg$cat3))
tg = cbind(tg, gene=!(tg$crp|tg$nbs|tg$te))

get_ovlp_len <- function(tq, tt) {
  gq = GRanges(seqnames=Rle(tq$chr), ranges=IRanges(tq$beg, end=tq$end))
  gqr = reduce(gq)
  gt = GRanges(seqnames=Rle(tt$chr), ranges=IRanges(tt$beg, end=tt$end))
  gtr = reduce(gt)
  
  qlen = sum(width(gqr))
  tlen = sum(width(gtr))
  olen = sum(width(intersect(gqr, gtr)))
  list('qlen'=qlen, 'tlen'=tlen, 'olen'=olen, 'pct'=olen/qlen)
}

fa = file.path(dir, "25.gal")
ta = read.table(fa, sep="\t", head=T, as.is=T)[,1:17]
fl = file.path(dir, "25.gall")
tl = read.table(fl, sep="\t", head=T, as.is=T)
fn = file.path(dir, "25.neti")
tn = read.table(fn, sep="\t", head=T, as.is=T)
sum(ta$id != tn$id)

tls = tl[tl$id %in% ta$id[tn$lev==1],]
get_ovlp_len(tg[tg$gene,], data.frame(chr=tls$tId, beg=tls$tBeg, end=tls$tEnd))

fd = file.path(dir, "27.tbl")
td = read.table(fd, sep="\t", head=F, as.is=T)
colnames(td) = c("chr", "beg", "end", "len")
gd = GRanges(seqnames=Rle(td$chr), ranges=IRanges(td$beg, end=td$end))

tn = tg[tg$crp,]
gn = GRanges(seqnames=Rle(tn$chr), ranges=IRanges(tn$beg, end=tn$end))
gdn = intersect(gd, gn)
tdn = data.frame(chr=as.character(seqnames(gdn)), beg=start(gdn), end=end(gdn), len=end(gdn)-start(gdn)+1)
tdn[order(tdn$len, decreasing=T),][1:20,]


chr1:4,254,121-4,261,352

chr8:3,905,901-3,907,283
chr5:41,783,926-41,786,169