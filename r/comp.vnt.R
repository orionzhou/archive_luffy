require(rtracklayer)
#require(PopGenome)
#require(pegas)
source("comp.fun.R")
source("Location.R")

tname = "HM101"
qnames = c(
  "HM058", "HM125", "HM056", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010"
)
dir = file.path(Sys.getenv("misc3"), "comp.vnt")

##### obtain regions covered by 9 (out of 12) accessions
tcfg = get_genome_cfg(tname)
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
  ti = ti[ti$lev <= 2,]
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

##### sliding window analysis
fw = file.path(Sys.getenv("misc3"), "compstat", "31.win.tbl")
tw = read.table(fw, header = T, sep = "\t", as.is = T)
gr = with(tw, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

fc = file.path(dir, "81.cvg.tbl")
tc = read.table(fc, header = F, sep = "\t", as.is = T)
grc = with(tc, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

lenc = intersect_basepair(gr, grc)

### calc nuc-div
cl = makeCluster(detectCores())
cluster_fun <- function() {}
clusterCall(cl, cluster_fun)

get_nuc_div <- function(rw, fs) {
  chr = rw['chr']; beg = as.numeric(rw['beg']); end = as.numeric(rw['end'])
  cmd = sprintf("tabix %s %s:%d-%d | cut -f6 | paste -sd+ | bc", fs, chr, beg, end)
  ret = as.numeric(system(cmd, intern = T))
  ifelse(is.numeric(ret), ret, NA)
}
fs = file.path(dir, '25.stat.tbl.gz')

ptm <- proc.time()
pis = parApply(cl, tw, 1, get_nuc_div, fs)
proc.time() - ptm

to = cbind(tw, lenc = lenc, pi = pis / lenc)
fo = file.path(Sys.getenv("misc3"), "compstat", "32.win.stat.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

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