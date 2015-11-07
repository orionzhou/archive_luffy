require(plyr)
require(dplyr)
require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(RColorBrewer)
source("Location.R")
source("comp.fun.R")

dirw = file.path(Sys.getenv("misc3"), "comp.sv.valid")

##### filter SVs (>= 50bp)
org = "HM034"
for (org in c("HM034", "HM056", "HM340")) {
dirg = file.path(Sys.getenv("genome"), org)

dirv = sprintf("%s/%s_HM101/31_sv", Sys.getenv("misc3"), org)
fv = file.path(dirv, "05.stb")
tv = read.table(fv, header = T, sep = "\t", as.is = T)
tv = tv[tv$tlen >= 50 | tv$qlen >= 50,]

fo = sprintf("%s/01_sv/%s.tbl", dirw, org)
write.table(tv, fo, sep = "\t", row.names = F, col.names = T, quote = F)
}

##### read in pacbio BAM for support (obsolete - see comp.sv.valid.py)
f_bam = sprintf("%s/pacbio/%s_%s/15.bam", Sys.getenv("misc3"), org, org)
bam = bamReader(f_bam, idx = T)

to = cbind(tv, nr1 = NA, nr2 = NA)
for (i in 1:100) {
#  i = 1
  chr = tv$qchr[i]; beg = tv$qbeg[i]; end = tv$qend[i]
  b1 = max(1, beg - 10)
  e1 = beg + 10
  b2 = max(1, end - 10)
  e2 = end + 10
  gr1 = GRanges(seqnames = chr, ranges = IRanges(b1, end = e1))
  gr2 = GRanges(seqnames = chr, ranges = IRanges(b2, end = e2))
  
  ds1 = read_bam(bam, gr1)
  ds2 = read_bam(bam, gr2)
  to$nr1[i] = nrow(ds1)
  to$nr2[i] = nrow(ds2)
}

sum(to$nr1 > 1 & to$nr2 > 1)
sum(to$nr1 <= 1 & to$nr2 <= 1)

####
for (org in c("HM034", "HM056", "HM340")) {
  fp = sprintf("%s/02_pacbio/%s.tbl", dirw, org)
  tp = read.table(fp, header = T, sep = "\t", as.is = T)
  cat(org, nrow(tp), sum(tp$n1 >= 5 & tp$n2 >= 5) / nrow(tp), "\n")
}
