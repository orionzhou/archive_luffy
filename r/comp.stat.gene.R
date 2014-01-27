library(Rsamtools)
source("comp.fun.R")

tname = "hm101"
qname = "hm056"
rname = "hm340"

t = read_genome_stat(tname)
q = read_genome_stat(qname)
r = read_genome_stat(rname)

vq = read_var_stat(qname)
vr = read_var_stat(rname)

cq = read_comp_stat(qname, tname)
cr = read_comp_stat(rname, tname)


tg = read.table(file.path(t$dir, "51.gtb")
  , sep="\t", as.is = T, header = T, quote="")[ , c(1, 3:6, 17)]
tg = cbind(tg, crp = grepl('CRP', tg$cat3), nbs = grepl('NBS', tg$cat3), 
  te = grepl('TE', tg$cat3))
tg = cbind(tg, gene=!(tg$crp|tg$nbs|tg$te))

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

# compare SV influences
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

