require(plyr)
require(dplyr)
require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(RColorBrewer)
source("Location.R")
source("comp.fun.R")

dirw = file.path(Sys.getenv('misc3'), 'comp.novseq')

##### compile contaminated scaffold IDs
qname = "HM340"
for (qname in qnames) {
cfg = cfgs[[qname]]
ccfg = ccfgs[[qname]]
tl = read.table(cfg$size, header = F, sep = "\t", as.is = T)
colnames(tl) = c("chr", "size")

fn = file.path(ccfg$cdir, "../41_novseq/12.foreign.bed")
tn = read.table(fn, header = F, sep = "\t", as.is = T)
colnames(tn) = c("chr", "beg", 'end', 'type')
tn$beg = tn$beg + 1
tn = cbind(tn, flen = tn$end - tn$beg + 1)

tm = ddply(tn, .(chr), summarise, flen = sum(flen))
tm = merge(tm, tl, by = 'chr')
tm = cbind(tm, fpct = tm$flen/tm$size)
fids = tm$chr[tm$fpct >= 0.5]
x1 = sum(tn$flen[tn$chr %in% fids]) / sum(tn$flen)
x2 = sum(tn$flen[!tn$chr %in% fids])

tg = read.table(cfg$gene, header = F, sep = "\t", as.is = T)
tg = tg[tg$V6 == 'mrna', c(1,5)]
colnames(tg) = c('chr', 'gid')
y = sum(tg$chr %in% fids)
cat(sprintf("%s: %d scfs rmvd [%d genes], %.03f nov-seq rmvd, %d bp left\n", qname, length(fids), y, x1, x2))

fo = file.path(ccfg$cdir, "../41_novseq/15.foreign.scf.txt")
write(fids, fo)
}

##### enrichment analysis
qname = "HM324"
cfg = cfgs[[qname]]
ccfg = ccfgs[[qname]]

tg = read.table(cfg$gene, sep = "\t", header = F, as.is = T)
colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
tg = tg[tg$type == 'cds',]
grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

fn = file.path(ccfg$cdir, "../41_novseq/21.bed")
tn = read.table(fn, sep = "\t", header = F, as.is = T)
grn = with(tn, GRanges(seqnames = V1, ranges = IRanges(V2+1, end = V3)))

olens = intersect_basepair(grg, grn)
ds = ddply(cbind(tg, olen = olens), .(fam), summarise, olen = sum(olen), alen = sum(end - beg + 1))
ds = cbind(ds, pct = ds$olen / ds$alen)
ds[order(ds$pct, decreasing = T),][1:20,]

