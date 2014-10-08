require(GenomicRanges)

org = "HM101"

chrs = sprintf("chr%s", 1:8)

dirr = file.path(Sys.getenv('genome'), org)
fa = file.path(dirr, '15.sizes')
fg = file.path(dirr, '16.gap.bed')
ta = read.table(fa, sep = '\t', header = F, as.is = T)
tg = read.table(fg, sep = '\t', header = F, as.is = T)
ta = ta[ta$V1 %in% chrs,]
tg = tg[tg$V1 %in% chrs,]
ga = GRanges(seqnames = ta$V1, ranges = IRanges(1, end = ta$V2))
gg = GRanges(seqnames = tg$V1, ranges = IRanges(tg$V2+1, end = tg$V3))
gr = setdiff(ga, gg)

fl = file.path(dirr, "51.tbl")
tl = read.table(fl, sep = '\t', header = F, as.is = T)
colnames(tl) = c("chr", "beg", "end", "srd", "id", "type", "cat")
tl = tl[tl$cat != 'TE',]

tc = tl[tl$type == 'cds',]
g_cds = GRanges(seqnames = tc$chr, ranges = IRanges(tc$beg, end = tc$end))
g_cds = reduce(g_cds)
g_cds = intersect(g_cds, gr)

tc = tl[tl$type == 'intron',]
g_ito = GRanges(seqnames = tc$chr, ranges = IRanges(tc$beg, end = tc$end))
g_ito = reduce(g_ito)
g_ito = setdiff(g_ito, g_cds)
g_ito = intersect(g_ito, gr)

tc = tl[tl$type == 'utr5' | tl$type == 'utr3',]
g_utr = GRanges(seqnames = tc$chr, ranges = IRanges(tc$beg, end = tc$end))
g_utr = reduce(g_utr)
g_utr = setdiff(g_utr, union(g_cds, g_ito))
g_utr = intersect(g_utr, gr)

g_ige = setdiff(gr, union(union(g_cds, g_ito), g_utr))


to = data.frame(chr = seqnames(go), beg = start(go) - 1, end = end(go))

types = c('CDS', 'Intron', 'UTR', 'Intergenic')
gos = list(g_cds, g_ito, g_utr, g_ige)

to = data.frame()
for (i in 1:length(types)) {
  type = types[i]
  go = gos[[i]]
  to = rbind(to, 
    data.frame(chr = seqnames(go), beg = start(go), end = end(go), type = type))
}

to = to[order(to$chr, to$beg, to$end), ]
fo = file.path(dirr, "51.merged.tbl")
write.table(to, fo, row.names = F, col.names = T, sep = "\t", quote = F)