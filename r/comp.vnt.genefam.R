require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(dplyr)
source("comp.fun.R")
source("Location.R")

dirw = file.path(Sys.getenv("misc3"), "comp.vnt")
diro = file.path(Sys.getenv("misc3"), "comp.stat")

tg = read.table(tcfg$gene, sep = "\t", header = F, as.is = T)
colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
tg = tg[tg$type == 'cds',]
grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

##### calculate ThetaPi for each gene
fa = file.path(dirw, "81.cvg.tbl")
ta = read.table(fa, header = F, sep = "\t", as.is = T)
gra = with(ta, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

fv = file.path(dirw, '25.stat.tbl')
tv = read.table(fv, header = T, sep = "\t", as.is = T)
gsnp = with(tv, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))

effmap = c(
  'downstream_gene_variant' = 'intergenic',
  'upstream_gene_variant' = 'intergenic',
  'intergenic_region' = 'intergenic',
  'missense_variant' = 'replacement',
  'synonymous_variant' = 'synonymous',
  'stop_retained_variant' = 'synonymous',
  'initiator_codon_variant' = 'synonymous',
  'intron_variant' = 'intron',
  'splice_region_variant' = 'intron',
  '5_prime_UTR_variant' = 'utr',
  '3_prime_UTR_variant' = 'utr',
  '5_prime_UTR_premature_start_codon_gain_variant' = 'utr',
  'stop_gained' = 'large_effect',
  'splice_acceptor_variant' = 'large_effect',
  'splice_donor_variant' = 'large_effect',
  'stop_lost' = 'large_effect',
  'start_lost' = 'large_effect'
)
ff = file.path(dirw, "23.snp.anno.tbl")
tf = read.table(ff, header = F, sep = "\t", as.is = T)
stopifnot(nrow(tv) == nrow(tf))
stopifnot(sum(!tf$V3 %in% names(effmap)) == 0)
tv = cbind(tv, eff = as.character(effmap[tf$V3]))


nds = intersect_score(grg, gsnp)
bps = intersect_basepair(grg, gra)

tg2 = cbind(tg, nd = nds, bp = bps)
gb = group_by(tg2, id)
dg = dplyr::summarise(gb, fam = fam[1], len = sum(end-beg+1), lenc = sum(bp), nd = sum(nd))

dg = dg[dg$lenc / dg$len >= 0.8,]
genes_covered = dg$id
dg = cbind(dg, pi = dg$nd / dg$lenc)
to = ddply(dg, .(fam), summarise, cnt = length(id), q25 = quantile(pi, 0.25), q50 = quantile(pi, 0.5), q75 = quantile(pi, 0.75))

fo = file.path(diro, "42.genefam.pi.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### characterize large-effect SNPs by gene fam
lefmap = c(
  'stop_gained' = 'premature stop',
  'splice_acceptor_variant' = 'splice site',
  'splice_donor_variant' = 'splice site',
  'stop_lost' = 'stop codon lost',
  'start_lost' = 'start codon lost'
)
idxs = tf$V3 %in% names(lefmap)
tvl = cbind(tv[idxs, c('chr','pos')], eff = as.character(lefmap[tf$V3[idxs]]))

dz = data.frame()
for (lef in lefmap) {
  tvls = tvl[tvl$eff == lef,]
  gvl = with(tvls, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))
  cnts = intersect_count(grg, gvl)
  dz = rbind(dz, cbind(tg, eff = lef, cnt = cnts))
}
gb = group_by(dz, id, eff)
dz2 = dplyr::summarize(gb, fam = fam[1], cnt = sum(cnt))
dg = dz2[dz2$id %in% genes_covered,]

to = ddply(dg, .(fam, eff), summarise, prop = sum(cnt > 0) / length(cnt))

fo = file.path(diro, "42.genefam.largeeff.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### Gene-Family Theta-Pi using reference mapping-based approach
dirw = file.path(Sys.getenv("misc3"), "hapmap/12_ncgr")
diro = file.path(Sys.getenv("misc3"), "comp.stat")

fa = file.path(dirw, "81.cvg.tbl")
ta = read.table(fa, header = F, sep = "\t", as.is = T)
gra = with(ta, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

fv = file.path(dirw, '45.stat.tbl')
tv = read.table(fv, header = T, sep = "\t", as.is = T)
gsnp = with(tv, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos), score = nucdiv))

nds = intersect_score(grg, gsnp)
bps = intersect_basepair(grg, gra)

tg2 = cbind(tg, nd = nds, bp = bps)
gb = group_by(tg2, id)
dg = dplyr::summarise(gb, fam = fam[1], len = sum(end-beg+1), lenc = sum(bp), nd = sum(nd))

dg = dg[dg$lenc / dg$len >= 0.8,]
genes_covered = dg$id
dg = cbind(dg, pi = dg$nd / dg$lenc)
to = ddply(dg, .(fam), summarise, cnt = length(id), q25 = quantile(pi, 0.25), q50 = quantile(pi, 0.5), q75 = quantile(pi, 0.75))

fo = file.path(diro, "42.genefam.pi.refmap.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

## restrict SNP calling to syntenic regions
diry = file.path(Sys.getenv("misc3"), "comp.vnt")
fa = file.path(diry, "81.cvg.tbl")
ta = read.table(fa, header = F, sep = "\t", as.is = T)
gry = with(ta, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

cat(sum(width(gra)), "\n")
cat(sum(width(gry)), "\n")
cat(sum(width(GenomicRanges::intersect(gra, gry))), "\n")

gra = GenomicRanges::intersect(gra, gry)
gsnp = GenomicRanges::intersect(gsnp, gry)

nds = intersect_score(grg, gsnp)
bps = intersect_basepair(grg, gra)

tg2 = cbind(tg, nd = nds, bp = bps)
gb = group_by(tg2, id)
dg = dplyr::summarise(gb, fam = fam[1], len = sum(end-beg+1), lenc = sum(bp), nd = sum(nd))

dg = dg[dg$lenc / dg$len >= 0.8,]
genes_covered = dg$id
dg = cbind(dg, pi = dg$nd / dg$lenc)
to = ddply(dg, .(fam), summarise, cnt = length(id), q25 = quantile(pi, 0.25), q50 = quantile(pi, 0.5), q75 = quantile(pi, 0.75))

fo = file.path(diro, "42.genefam.pi.refmap.syn.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
