require(rtracklayer)
require(GenomicRanges)
#require(PopGenome)
#require(pegas)
source("comp.fun.R")
source("Location.R")

tname = "HM101"
qnames = get_orgs('ingroup')

dirw = file.path(Sys.getenv("misc3"), "comp.vnt")
diro = file.path(Sys.getenv("misc3"), "comp.stat")

##### obtain regions covered by 9 (out of 12) accessions
tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

gr = setdiff(grt, grp)
gra = GRanges()
for (qname in qnames) {
  diri = sprintf("%s/%s_%s", Sys.getenv("misc3"), qname, tname)
  
  fi = file.path(diri, '23_blat/31.9/gax')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:3,10)]
  colnames(ti) = c('chr', 'beg', 'end', 'lev')
  ti = ti[ti$lev <= 1,]
  grg = with(ti, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
  gra = c(gra, grg)
}
gr_cvg = as(coverage(gra), "GRanges")
for (i in 1:length(qnames)) {
  cat(sprintf("%2d: %d\n", i, sum(width(gr_cvg[mcols(gr_cvg)$score >= i, ]))))
}
grc = reduce(gr_cvg[mcols(gr_cvg)$score >= 9, ])

tcvg = data.frame(chr = seqnames(grc), beg = start(grc), end = end(grc))
fo = file.path(dir, "81.cvg.tbl")
write.table(tcvg, fo, sep = "\t", row.names = F, col.names = F, quote = F)

##### creating sliding window statistics
fw = file.path(Sys.getenv("misc3"), "compstat", "31.win.tbl")
tw = read.table(fw, header = T, sep = "\t", as.is = T)
gr = with(tw, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

fc = file.path(dir, "81.cvg.tbl")
tc = read.table(fc, header = F, sep = "\t", as.is = T)
grc = with(tc, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

lenc = intersect_basepair(gr, grc)

### calc nuc-div for SNPs
fn = file.path(Sys.getenv("misc3"), "comp.vnt", "25.stat.tbl")
tn = read.table(fn, header = T, sep = "\t", as.is = T)
grn = with(tn, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))
nucdivs = intersect_score(gr, grn)
pi_snp = nucdivs / lenc

# (obsolete!) calc nuc-div for SNPs - the hard way
cl = makeCluster(detectCores())
cluster_fun <- function() {}
clusterCall(cl, cluster_fun)

get_nuc_div <- function(rw, fs) {
  chr = rw['chr']; beg = as.numeric(rw['beg']); end = as.numeric(rw['end'])
  cmd = sprintf("tabix %s %s:%d-%d | cut -f6 | paste -sd+ | bc", fs, chr, beg, end)
  ret = as.numeric(system(cmd, intern = T))
  ifelse(is.numeric(ret), ret, NA)
}
fs = file.path(dirw, '25.stat.tbl.gz')

ptm <- proc.time()
pis = parApply(cl, tw, 1, get_nuc_div, fs)
proc.time() - ptm

### calc nuc-div for indels/sv
fv = file.path(Sys.getenv("misc3"), "comp.vnt", "52.stat.tbl")
tv = read.table(fv, header = T, sep = "\t", as.is = T)
tv = cbind(tv, size = (tv$rsize + tv$asize - 2))
tvs = tv[tv$size < 50,]
tvl = tv[tv$size >= 50,]
grvs = with(tvs, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))
grvl = with(tvl, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))

nucdivs = intersect_score(gr, grvs)
pi_indel = nucdivs / tw$len_ng
nucdivs = intersect_score(gr, grvl)
pi_sv = nucdivs / tw$len_ng

# output
pi_snp[is.infinite(pi_snp)] = NA
pi_indel[is.infinite(pi_indel)] = NA
pi_sv[is.infinite(pi_sv)] = NA
to = cbind(tw, lenc = lenc, pi_snp = pi_snp, pi_indel = pi_indel, pi_sv = pi_sv)
fo = file.path(diro, "32.win.stat.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

##### misc tests of packages
#popgenome
dv = file.path(dir, "xxx")
x1 = readData(dv, include.unknown = T, format = "VCF")
get.sum.data(x1)

x1 = detail.stats(x1)
mafs = x1@region.stats@minor.allele.freqs[[1]]

x1 = F_ST.stats(x1)
x1 = diversity.stats(x1)
x1@nuc.diversity.within

# variantannotation + snpStats
require(VariantAnnotation)
require(snpStats)

fv = file.path(dir, "xxx/x.vcf")
vcf <- readVcf(fv, "mt40")
res = genotypeToSnpMatrix(vcf)
cat(nrow(res$map), "sites", sum(res$map$ignore), "poly-allelic\n", sep = " ")
snpsum = col.summary(res$genotypes)[!res$map$ignore,]
sum( (2*snpsum$Calls)/(2*snpsum$Calls-1) * 2 * snpsum$MAF * (1-snpsum$MAF) )

##### genome-wide population genetics statistics
fa = file.path(dirw, "81.cvg.tbl")
ta = read.table(fa, header = F, sep = "\t", as.is = T)
gra = with(ta, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

f1 = file.path(Sys.getenv("genome"), "HM101", "51.tbl")
t1 = read.table(f1, header = F, sep = "\t", as.is = T)
colnames(t1) = c("chr", "beg", "end", "srd", "id", "type", "fam")

tt = t1[t1$type == 'cds',]
gcds = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gcds = GenomicRanges::intersect(gra, reduce(gcds))

tt = t1[t1$type == 'intron',]
gito = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gito = GenomicRanges::intersect(GenomicRanges::setdiff(gra, gcds), reduce(gito))

tt = t1[t1$type == 'utr5',]
gut5 = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gut5 = GenomicRanges::intersect(GenomicRanges::setdiff(gra, c(gcds, gito)), reduce(gut5))

tt = t1[t1$type == 'utr3',]
gut3 = with(tt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
gut3 = GenomicRanges::intersect(GenomicRanges::setdiff(gra, c(gcds, gito, gut5)), reduce(gut3))

gitr = GenomicRanges::setdiff(gra, reduce(c(gcds, gito, gut5, gut3)))

f2 = file.path(Sys.getenv("genome"), "HM101", "51.syn.tbl")
t2 = read.table(f2, header = F, sep = "\t", as.is = T)
gsyn = with(t2, GRanges(seqnames = V1, ranges = IRanges(V2, end = V2)))
gsyn = GenomicRanges::intersect(gra, reduce(gsyn))
grpl = GenomicRanges::setdiff(gcds, gsyn)

grl = GRangesList(cds = gcds, synonymous = gsyn, replacement = grpl, intron = gito, utr5 = gut5, utr3 = gut3, intergenic = gitr)

fn = file.path(dirw, '25.stat.tbl')
tn = read.table(fn, header = T, sep = "\t", as.is = T)
gsnp = with(tn, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))

  pi = sum(tn$nucdiv) / sum(width(gra))
  tw = table(tn$nsam)
  ans = sapply(as.numeric(names(tw)), get_har <- function(x) sum(1 / 1:(x-1)))
  thetaw = sum(as.numeric(tw) / ans) / sum(width(gra))

stats = list()
stats[["Total"]] = matrix(c(
    format(sum(width(gra)), big.mark = ","), 
    '-', 
    format(nrow(tn), big.mark = ","),
    sprintf("%.04f", pi), 
    sprintf("%.04f", thetaw)), nrow = 1, 
    dimnames = list(NULL, c("Covered bases", "Percent",
    "Polymorphic sites", "Pi", "ThetaW")))
for (i in 1:length(grl)) {
  ma = as.matrix(findOverlaps(gsnp, grl[i]))
  idxs = ma[,1]
  tns = tn[idxs,]
  
  bp = sum(width(grl[[i]]))
  pct = bp / sum(width(gra))
  n_snp = nrow(tns)
  pi = sum(tns$nucdiv) / bp
  
  tw = table(tns$nsam)
  ans = sapply(as.numeric(names(tw)), get_har <- function(x) sum(1 / 1:(x-1)))
  thetaw = sum(as.numeric(tw) / ans) / bp

  stat = c(
    format(bp, big.mark = ","), 
    sprintf("%.02f", pct), 
    format(n_snp, big.mark = ","),
    sprintf("%.04f", pi), 
    sprintf("%.04f", thetaw))
  stats[[names(grl)[i]]] = matrix(stat, nrow = 1, 
    dimnames = list(NULL, c("Covered bases", "Percent",
    "Polymorphic sites", "Pi", "ThetaW")))
}
ds = do.call(rbind.data.frame, stats)
fo = file.path(diro, "13_diversity.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)

