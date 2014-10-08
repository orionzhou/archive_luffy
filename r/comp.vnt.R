require(GenomicRanges)
require(plyr)
require(rtracklayer)
require(Rsamtools)
require(ggplot2)
require(grid)
#require(VennDiagram)

orgs = c(
  "HM058", "HM056", "HM125", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM010", "HM018", "HM022", "HM324", "HM340"
)
#orgs = toupper(orgs)

##### Mapping-based approach
chrs = sprintf("chr%s", 1:8)

dirr = file.path(Sys.getenv('genome'), "HM101")
fa = file.path(dirr, '15.sizes')
fg = file.path(dirr, '16.gap.bed')
ta = read.table(fa, sep = '\t', header = F, as.is = T)
tg = read.table(fg, sep = '\t', header = F, as.is = T)
ta = ta[ta$V1 %in% chrs,]
tg = tg[tg$V1 %in% chrs,]
ga = GRanges(seqnames = ta$V1, ranges = IRanges(1, end = ta$V2))
gg = GRanges(seqnames = tg$V1, ranges = IRanges(tg$V2+1, end = tg$V3))
gr = setdiff(ga, gg)

ff = file.path(dirr, "51.tbl")
tf = read.table(ff, sep = '\t', header = F, as.is = T)
colnames(tf) = c("chr", "beg", "end", "srd", "id", "type", "cat")

fc = file.path(dirr, "51.merged.tbl")
tc = read.table(fc, sep = '\t', header = T, as.is = T)
gc = GRanges(seqnames = tc$chr, ranges = IRanges(tc$beg, end = tc$end))

dirm = file.path(Sys.getenv("misc3"), "hapmap", "12_ncgr")

### generate coverage BED file (takes long time - run with caution)
for (org in orgs) {
fcov = sprintf("%s/35_cov/%s.bw", dirm, org)

bw = import(fcov, which = gr, asRangedData = F)
chrs = seqnames(bw)
poss = start(bw)
idxs = bw$score >= 1
rm(bw)
gm = GRanges(seqnames = chrs[idxs], ranges = IRanges(poss[idxs], 
  end = poss[idxs]))
gm = reduce(gm)
rm(chrs, poss, idxs)

tm = data.frame(chr = seqnames(gm), beg = start(gm) - 1, end = end(gm))
fo = sprintf("%s/38_covered/%s.bed", dirm, org)
write.table(tm, fo, row.names = F, col.names = F, sep = "\t", quote = F)
}

### SNP density
tr = data.frame()

for (org in orgs) {
#org = "HM004"

fs = file.path(dirm, "42_vnt", org, "snp")
ts = read.table(fs, sep = "\t", header = F, as.is = T)
colnames(ts) = c("chr", "pos", "ref", "alt")
gs = GRanges(seqnames = ts$chr, ranges = IRanges(ts$pos, end = ts$pos))

fv = sprintf("%s/38_covered/%s.bed", dirm, org)
tv = read.table(fv, sep = "\t", header = F, as.is = T)
gv = GRanges(seqnames = tv$V1, ranges = IRanges(tv$V2+1, end = tv$V3))

for (type in unique(tc$type)) {
  gcs = gc[tc$type == type]
  gsi = intersect(gs, gcs)
  nsnp = sum(width(gsi))
  
  gvi = intersect(gv, gcs)
  bcov = sum(width(gvi))
  trs = data.frame(org = org, type = type, bcov = bcov, nsnp = nsnp)
  tr = rbind(tr, trs)
}
}

tr = cbind(tr, subrate = tr$nsnp / tr$bcov)

##### Assembly-based approach
### SNP density
ta = data.frame()

for (org in orgs) {
#org = "HM034"
dira = sprintf("%s/%s_HM101/23_blat/31.9", Sys.getenv('misc3'), org)

fs = file.path(dira, "snp")
ts = read.table(fs, sep = "\t", header = F, as.is = T)
colnames(ts) = c("chr", "pos", "ref", "alt", "qid", "qpos", "cid", "lev")
ts = ts[ts$lev <= 2,]
gs = GRanges(seqnames = ts$chr, ranges = IRanges(ts$pos, end = ts$pos))

fv = file.path(dira, "gax")
tv = read.table(fv, sep = "\t", header = F, as.is = T)
colnames(tv) = c("chr", "beg", "end", "srd", "qid", "qbeg", "qend", "qsrd",
  "cid", "lev")
tv = tv[tv$lev <= 2,]
gv = GRanges(seqnames = tv$chr, ranges = IRanges(tv$beg, end = tv$end))

for (type in unique(tc$type)) {
  gcs = gc[tc$type == type]
  gsi = intersect(gs, gcs)
  nsnp = sum(width(gsi))
  
  gvi = intersect(gv, gcs)
  bcov = sum(width(gvi))
  tas = data.frame(org = org, type = type, bcov = bcov, nsnp = nsnp)
  ta = rbind(ta, tas)
}
}

ta = cbind(ta, subrate = ta$nsnp / ta$bcov)


##### compare SNPs called from two sources
to = rbind(cbind(method = 'mapping', tr), cbind(method = 'assembly', ta))
fo = file.path(Sys.getenv("misc3"), "compstat", "snp.pdf")

to$type = factor(to$type, levels = c('CDS', 'Intron', 'UTR', 'Intergenic'))
to$method = factor(to$method, levels = c("mapping", "assembly"))
p = ggplot(to) +
  geom_bar(mapping = aes(x = type, y = subrate, fill = method), 
    stat = 'identity', position = 'dodge', geom_params=list(width = 0.5)) +
  scale_fill_brewer(palette='Set2', name = '',
#    breaks = c('assembly', 'mapping'), 
    labels = c("Mapping-based", "Assembly-based")) +
  scale_x_discrete(name = '') +
  scale_y_continuous(name = 'SNP density (/bp)') +
  facet_wrap( ~ org, nrow = 5) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top", legend.key.size = unit(0.5, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(0,1,1,0), "lines")) +
  theme(axis.text.x = element_text(size = 8, colour = "grey", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "blue", angle = 0)) +
  theme(panel.margin = unit(0.6, 'lines'))
ggsave(p, filename = fo, width = 6, height = 6)


### venn diagram
chrs = seq(8)
names(chrs) = sprintf("chr%s", seq(8))

dir = sprintf("%s/%s_%s/23_blat", DIR_Misc3, qorg, torg)
fd = file.path(dir, "27.snp")
td = read.table(fd, sep='\t', header=T, as.is=T)[,3:4]
colnames(td) = c('chr', 'pos')
td.1 = td[td$chr %in% names(chrs),]
possd = chrs[td.1$chr]*1000000000+td.1$pos

fb = sprintf('/home/youngn/zhoup/Data/misc3/hapmap_mt40/30_vnt/%s.snp', qorg)
tb = read.table(fb, sep='\t', header=T, as.is=T)[,1:2]
colnames(tb) = c('chr', 'pos')
tb.1 = tb[tb$chr %in% names(chrs),]
possb = chrs[tb.1$chr]*1000000000+tb.1$pos

area1 = length(possd)
area2 = length(possb)
areac = sum(possd %in% possb)
venn.plot <- draw.pairwise.venn( area1, area2, areac, 
category = c("de novo", "read mapping"),
fill = c("blue", "red"), lty = "blank", cex = 2, 
cat.cex = 2, cat.pos = c(285, 105), cat.dist = 0.09, cat.just = list(c(-1, -1), c(1, 1)),
ext.pos = 30, ext.dist = -0.05, ext.length = 0.85, ext.line.lwd = 2, ext.line.lty = "dashed")
tiff(sprintf('/home/youngn/zhoup/Data/misc3/%s_%s/snpcmp.tiff', qorg, torg))
grid.draw(venn.plot)
dev.off()


