require(GenomicRanges)
require(plyr)
require(rtracklayer)
#require(VennDiagram)

org = "HM034"

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

### generate coverage BED file (takes long time - run with caution)
orgs = c(
  "hm004", "hm010", "hm018", "hm022", "hm034", 
  "hm050", "hm056", "hm058", "hm060", "hm095", 
  "hm125", "hm129", "hm185", "hm324", "hm340")
orgs = c(
  "hm050", "hm056", "hm058", "hm060", "hm095", 
  "hm125", "hm129", "hm185", "hm324", "hm340")
orgs = toupper(orgs)

for (org in orgs) {
dirm = file.path(Sys.getenv("misc3"), "hapmap", "12_ncgr")
fcov = sprintf("%s/35_cov/%s.bw", dirm, org)

bw = import(fcov, which = gr, asRangedData = F)
chrs = seqnames(bw)
poss = start(bw)
idxs = bw$score >= 1
rm(bw)
gm = GRanges(seqnames = chrs[idxs], ranges = IRanges(poss[idxs], end = poss[idxs]))
gm = reduce(gm)
rm(chrs, poss, idxs)

tm = data.frame(chr = seqnames(gm), beg = start(gm) - 1, end = end(gm))
fo = sprintf("%s/38_covered/%s.bed", dirm, org)
write.table(tm, fo, row.names = F, col.names = F, sep = "\t", quote = F)
}

# SNP density
fw = sprintf("%s/%s_HM101/23_blat/31.9/snp", Sys.getenv("misc3"), org)
tw = read.table(fw, sep='\t', header=T, as.is=T)[,1:17]
tw = cbind(tw, alnlen=tw$match+tw$misMatch)
tw = tw[tw$alnlen > 0,]

sum(tw$alnlen)
sum(tw$misMatch)
sum(tw$misMatch) / sum(tw$alnlen)

intvs = c(0,100,1000,10000,100000,Inf)
labels = c("1-100", "100-1000", "1k-10k", "10k-100k", "100k+")

tw.2 = cbind(tw, sub=tw$misMatch/tw$alnlen, intv=cut(tw$alnlen, breaks=intvs, labels=labels))
tws = ddply(tw.2, .(intv), summarise, cnt=length(sub), misMatch=sum(misMatch), alnlen=sum(alnlen), sub.mean=sum(misMatch)/(sum(match)+sum(misMatch)), sub.median=median(sub))

# compare SNPs called from two sources
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


