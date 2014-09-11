require(rtracklayer)
source("comp.fun.R")

tname = "hm101"
qname1 = "hm056"
qname2 = "hm034"
qname3 = "hm340"
qname = qname1

t = read_genome_stat(tname)
q = read_genome_stat(qname)
cq = read_comp_stat(qname, tname)
#vq = read_var_stat(qname)


tg = read.table(file.path(t$dir, "51.gtb")
  , sep="\t", as.is = T, header = T, quote="")[ , c(1, 3:6, 17)]
tg = cbind(tg, crp = grepl('CRP', tg$cat3), nbs = grepl('NBS', tg$cat3), 
  te = grepl('TE', tg$cat3))
tg = cbind(tg, gene=!(tg$crp | tg$nbs | tg$te))

tgg = tg[tg$gene == T, ]
tgt = tg[tg$te == T, ]
tgc = tg[tg$crp == T, ]
tgn = tg[tg$nbs == T, ]
tg = rbind(tgg[sample(nrow(tgg), 1000), ], tgt[sample(nrow(tgt), 500), ],
  tgc, tgn)

gg = GRanges(seqnames = Rle(tg$chr),
  ranges = IRanges(tg$beg, end = tg$end))
glg = with(tg, makeGRangesListFromFeatureFragments(
  seqnames = chr, fragmentStarts = sprintf("%d,", beg), 
  fragmentWidths = sprintf("%d,", end - beg + 1), strand = srd))

grg = GRanges(seqnames = tgg$chr, ranges = IRanges(tgg$beg, end = tgg$end))
grt = GRanges(seqnames = tgt$chr, ranges = IRanges(tgt$beg, end = tgt$end))
grc = GRanges(seqnames = tgc$chr, ranges = IRanges(tgc$beg, end = tgc$end))
grn = GRanges(seqnames = tgn$chr, ranges = IRanges(tgn$beg, end = tgn$end))

grlg = with(tgg, makeGRangesListFromFeatureFragments(
  seqnames = chr, fragmentStarts = sprintf("%d,", beg), 
  fragmentWidths = sprintf("%d,", end - beg + 1), strand = srd))
grlt = with(tgt, makeGRangesListFromFeatureFragments(
  seqnames = chr, fragmentStarts = sprintf("%d,", beg), 
  fragmentWidths = sprintf("%d,", end - beg + 1), strand = srd))
grlc = with(tgc, makeGRangesListFromFeatureFragments(
  seqnames = chr, fragmentStarts = sprintf("%d,", beg), 
  fragmentWidths = sprintf("%d,", end - beg + 1), strand = srd))
grln = with(tgn, makeGRangesListFromFeatureFragments(
  seqnames = chr, fragmentStarts = sprintf("%d,", beg), 
  fragmentWidths = sprintf("%d,", end - beg + 1), strand = srd))

# assess conserved proportion of genic regions
fgax = file.path(cq$dir, "31.9/gax")
tgax = read.table(fgax, sep = "\t", header = F, as.is = T)
colnames(tgax) = c("tId", "tBeg", "tEnd", "tSrd", "id", "qId", "qBeg", "qEnd",
  "qSrd")
gr = GRanges(seqnames = tgax$tId, ranges = IRanges(tgax$tBeg, end = tgax$tEnd))
grl = with(tgax, makeGRangesListFromFeatureFragments(
  seqnames = tId, fragmentStarts = sprintf("%d,", tBeg), 
  fragmentWidths = sprintf("%d,", tEnd - tBeg + 1), strand = tSrd))

get_ovlp_len <- function(dfi, gr) {
  c('leno' = sum(width(gr[dfi$tidx])))
}
get_ovlp_90_idxs <- function(tgr, tgrl, qgr, qgrl) {
  gov = intersect(qgr, tgr)
  ma = as.matrix(findOverlaps(gov, qgrl))
  dfi = data.frame(idx = ma[,2], tidx = ma[,1])
  dfl = ddply(dfi, .(idx), get_ovlp_len, gov)
  df0 = data.frame(idx = 1:length(qgr), len = width(qgr))
  df = merge(df0, dfl, by = 'idx', all = T)
  df$leno[is.na(df$leno)] <- 0
  df = cbind(df, pct = df$leno / df$len)
  
  n80 = sum(df$pct >= 0.8)
  n90 = sum(df$pct >= 0.9)
  cat(80, n80, n80/nrow(df), "\n", sep = "\t")
  cat(90, n90, n90/nrow(df), "\n", sep = "\t")
  df
}
dog = get_ovlp_90_idxs(gr, grl, grg, grlg)
dot = get_ovlp_90_idxs(gr, grl, grt, grlt)
don = get_ovlp_90_idxs(gr, grl, grn, grln)
doc = get_ovlp_90_idxs(gr, grl, grc, grlc)

#

gr_idm_l = GRanges(seqnames = tm$tId, 
  ranges = IRanges(tm$tBeg, end = tm$tEnd))

get_ovlp_idxs <- function(grl, gr) {
  ma = as.matrix(findOverlaps(gr, grl))
  didx1 = data.frame(sidx=ma[,2], qidx=ma[,1])
  idxs_rm = unique(didx1$sidx)
  cat(length(idxs_rm), length(idxs_rm) / length(grl), "\n", sep="\t")
  idxs_rm
}
idxg = get_ovlp_idxs(grlg, gr_idm_l)
idxt = get_ovlp_idxs(grlt, gr_idm_l)
idxc = get_ovlp_idxs(grlc, gr_idm_l)
idxn = get_ovlp_idxs(grln, gr_idm_l)

tmg = tgg[idxg,]
tmt = tgt[idxt,]
tmc = tgc[idxc,]
tmn = tgn[idxn,]

tgm = rbind(tmg, tmt, tmc, tmn)
lent = tgm$end - tgm$beg + 1
tgm = cbind(tgm, loc = sprintf("%s:%d-%d", tgm$chr, tgm$beg - lent/2, 
  tgm$end + lent/2))
write.table(tgm, file.path(dir, "81.gene.tbl"), sep = "\t",
  quote = F, row.names = F, col.names = T)


# compare SNP/indel density
ds = vq

snp = ds$snp
gs = GRanges(seqnames = Rle(snp$chr), 
    ranges = IRanges(snp$pos, end = snp$pos))

ma = as.matrix(findOverlaps(gs, glg))
didx1 = data.frame(sidx=ma[,2], qidx=ma[,1])
didx <- ddply(didx1, .(sidx), nrow)
colnames(didx)[2] = "nsnp"

tgn = merge(cbind(tg, sidx = 1 : nrow(tg)), didx, by = 'sidx', all = T)
tgn$nsnp[is.na(tgn$nsnp)] = 0
tgn = cbind(tgn, snpdensity = tgn$nsnp / (tgn$end - tgn$beg + 1))

sum_by_col <- function(cname, cname.stat, df) {
    da = df[df[ ,cname] == T, cname.stat]
    list(mean = mean(da), median = median(da), sd = sd(da))
}
sapply(c('gene', 'te', 'crp', 'nbs'), sum_by_col, 'snpdensity', tgn)

# compare INS/DEL density
tid = read.table(file.path(cq$dir, "23_blat/27.vnt/05.tbl"), 
  header = F, sep = "\t", as.is = T)
colnames(tid) = c('tid', 'tbeg', 'tend', 'id', 'qid', 'qbeg', 'qend', 
  'qlen', 'tlen')
gri = GRanges(seqnames = tid$tid, ranges = IRanges(tid$tbeg, end = tid$tend))

# compare coverage
ds = vq

get_pct_covg <- function(x, bamFile) {
  chr = as.character(x['chr'])
  beg = as.integer(x['beg'])
  end = as.integer(x['end'])
  len = end - beg + 1
  param <- ScanBamParam(what = c('pos', 'qwidth'),
    which = GRanges(chr, IRanges(beg, end)),
    flag = scanBamFlag(isUnmappedQuery = F, isDuplicate = F))
  y <- scanBam(bamFile, param=param)[[1]]
  covg = coverage(IRanges(y[["pos"]], width=y[["qwidth"]]))
  if(length(covg) < end) {
    covg = c(covg, Rle(0, end-length(covg)))
  }
  sum(covg[beg:end] >= 3) / len
}

ptm <- proc.time()
covg = apply(tg[1:1000, ], 1, get_pct_covg, ds$fbam)
proc.time() - ptm


# compare blat coverage
ds = cq
tw = ds$tw
g_cov <- reduce(GRanges(seqnames = Rle(tw$tId),
  ranges = IRanges(tw$tBeg, end = tw$tEnd)))

ma = as.matrix(findOverlaps(g_cov, glg))
didx = data.frame(sidx=ma[,2], qidx=ma[,1])

get_len_ovlp <- function(df, glg, gr) {
  sidx = df$sidx[1]
#  gg = glg[[sidx]]
  grs = gr[df$qidx]
  data.frame(sidx = sidx, len_ovlp = sum(width(grs)))
}
tcov <- ddply(didx, .(sidx), get_len_ovlp, glg, g_cov)

tgn = merge(cbind(tg, sidx = 1 : nrow(tg)), tcov, by = 'sidx', all = T)
tgn$len_ovlp[is.na(tgn$len_ovlp)] = 0
tgn = cbind(tgn, pct_cov = tgn$len_ovlp / (tgn$end - tgn$beg + 1))

sapply(c('gene', 'te', 'crp', 'nbs'), sum_by_col, 'pct_cov', tgn)

# compare SV impacts
indel = ds$indel
g_sv <- reduce(GRanges(seqnames = Rle(indel$tid),
  ranges = IRanges(indel$tbeg+1, end = indel$tend)))
g_sv_g <- intersect(gg, g_sv)

ma = as.matrix(findOverlaps(g_sv_g, glg))
didx = data.frame(sidx=ma[,2], qidx=ma[,1])

tsv <- ddply(didx, .(sidx), get_len_ovlp, glg, g_sv)

tgn = merge(cbind(tg, sidx = 1 : nrow(tg)), tsv, by = 'sidx', all = T)
tgn$len_ovlp[is.na(tgn$len_ovlp)] = 0
tgn = cbind(tgn, pct_sv = tgn$len_ovlp / (tgn$end - tgn$beg + 1))

sapply(c('gene', 'te', 'crp', 'nbs'), sum_by_col, 'pct_sv', tgn)

