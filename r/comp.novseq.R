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
diro = file.path(Sys.getenv('misc3'), 'comp.stat')

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
do = data.frame()
for (qname in qnames_all[1:12]) {
cfg = cfgs[[qname]]

tg = read.table(cfg$gene, sep = "\t", header = F, as.is = T)
colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "fam")
tg = tg[tg$type == 'cds',]
grg = with(tg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

fn = file.path(cfg$cdir, "../41_novseq/21.bed")
tn = read.table(fn, sep = "\t", header = F, as.is = T)
grn = with(tn, GRanges(seqnames = V1, ranges = IRanges(V2+1, end = V3)))

olens = intersect_basepair(grg, grn)
#ds = ddply(cbind(tg, olen = olens), .(fam), summarise, olen = sum(olen), alen = sum(end - beg + 1))
#ds = cbind(ds, prop = ds$olen / ds$alen, org = qname)
gb = group_by(cbind(tg, olen = olens), id)
ds1 = summarise(gb, len = sum(end-beg+1), olen = sum(olen), pct = olen/len, fam = fam[1])
ds = ddply(ds1, .(fam), summarise, cnt = length(pct), prop = sum(pct>=0.5)/cnt)

do = rbind(do, cbind(ds, org = qname))
}
to = ddply(do, .(fam), summarise, q25 = quantile(prop, 0.25), q50 = quantile(prop, 0.5), q75 = quantile(prop, 0.75))

fo = file.path(diro, "42.genefam.novseq.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)


fi = file.path(diro, "42.genefam.novseq.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)

ffam = file.path(diro, "41.gene.fams.tbl")
fams = read.table(ffam, header = F, sep = "\t", as.is = T)[,1]

tis = ti[ti$fam %in% fams,]
tis = tis[order(tis$q50, tis$q25, tis$q75, decreasing = T),]
tis$fam = factor(tis$fam, levels = tis$fam)
p3 = ggplot(tis) +
  geom_crossbar(aes(x = fam, y = q50, ymin = q25, ymax = q75),
    stat = 'identity', position = 'dodge', geom_params = list(width = 0.7, size = 0.3)) + 
  coord_flip() +
  scale_x_discrete(name = '', expand = c(0.01, 0.01)) +
  scale_y_continuous(name = 'Proportion novel (absent in HM101)', expand = c(0.002, 0.002)) +
  theme_bw() +
#  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = file.path(diro, "18_novseq_genefam.pdf")
ggsave(p3, filename = fp, width = 5, height = 4)
