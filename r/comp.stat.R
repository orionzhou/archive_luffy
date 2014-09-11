require(rtracklayer)
require(plyr)
require(GenomicRanges)
source("comp.fun.R")

tname = "hm101"
qname = "hm340"

t = read_genome_stat(tname)
q = read_genome_stat(qname)
#cq = read_comp_stat(qname, tname)
#vq = read_var_stat(qname)

# basic assembly stats
library(Biostrings)
source(paste("http://faculty.ucr.edu/~tgirke/",
    "Documents/R_BioCond/My_R_Scripts/contigStats.R", sep=''))

scf = sprintf("/home/youngn/zhoup/Data/genome/%s/scf.fas", qname)
assembly <- readDNAStringSet(scf, "fasta")
N <- list(acc = width(assembly))
reflength <- sapply(N, sum)
stats <- contigStats(N = N, reflength = reflength, style = "data")
stats[["Contig_Stats"]]

ctg = sprintf("/home/youngn/zhoup/Data/genome/%s/ctg.fas", org)
assembly <- readDNAStringSet(ctg, "fasta")
N <- list(acc = width(assembly))
reflength <- sapply(N, sum)
stats <- contigStats(N = N, reflength = reflength, style = "data")
stats[["Contig_Stats"]]

# plot scaffold size distribution
tmp = cut(t_len$length / 1000, breaks = c(0,1,5,10,50,100,500,1000,5000))
p = ggplot(data.frame(size = tmp)) +
  geom_bar(aes(x = factor(size)), width = 0.7) + 
  scale_x_discrete(name = "Scaffold Size (kb)") + 
  scale_y_continuous(name = "") +
  theme(axis.text.x = element_text(angle = 15, size = 8))
ggsave(file.path(q$dir, "figs/01_scaffold_size.png"), p, width=5, height=4)


##### comp stats
qname = "hm034"
tname = "hm101"
dir = sprintf("/home/youngn/zhoup/Data/misc3/%s_%s/23_blat", 
  toupper(qname), toupper(tname))
fi = file.path(dir, '31.5.gal')
ti = read.table(fi, header = T, sep = "\t", as.is = T)[, c(1:19)]
ddply(ti, .(lev), summarise, ali=sum(ali))
fi = file.path(dir, '31.9.gal')
ti = read.table(fi, header = T, sep = "\t", as.is = T)[, c(1:19)]
ddply(ti, .(lev), summarise, ali=sum(ali))

fx = file.path(dir, '31.9/gax')
tx = read.table(fx, header = F, sep = "\t", as.is = T)
colnames(tx) = c("tid", 'tb', 'te', 'tsrd', 'qid', 'qb', 'qe', 'qsrd', 'cid', 'lev')
grq = GRanges(seqnames = tx$qid, ranges = IRanges(tx$qb, end = tx$qe))
sum(width(grq))
sum(width(reduce(grq)))



fgt = file.path(dir, '31.9/gax')
fgq = file.path(dir, '41.9/gax')
tt = read.table(fgt, header = F, sep = "\t", as.is = T)
colnames(tt) = c('tid', 'tbeg', 'tend', 'tsrd', 'qid', 'qbeg', 'qend', 'qsrd', 'cid', 'lev')
tq = read.table(fgq, header = F, sep = "\t", as.is = T)
colnames(tq) = c('qid', 'qbeg', 'qend', 'qsrd', 'tid', 'tbeg', 'tend', 'tsrd', 'cid', 'lev')

grt1 = GRanges(seqnames = tt$tid, ranges = IRanges(tt$tbeg, end = tt$tend))
grt2 = GRanges(seqnames = tq$tid, ranges = IRanges(tq$tbeg, end = tq$tend))

grq1 = GRanges(seqnames = tt$qid, ranges = IRanges(tt$qbeg, end = tt$qend))
grq2 = GRanges(seqnames = tq$qid, ranges = IRanges(tq$qbeg, end = tq$qend))
sum(width(grq1))
sum(width(grq2))
sum(width(intersect(grq1, grq2)))


fd = file.path(dir, 't.t.4.del.bed')
td = read.table(fd, header = F, sep = "\t", as.is = T)
colnames(td) = c('chr', 'beg', 'end')
td$beg = td$beg + 1
td = cbind(td, len = td$end - td$beg + 1, 
  str = sprintf("%s:%d-%d", td$chr, td$beg, td$end))
td[order(td$len, decreasing = T), ][1:50,]


# plot global pairwise comparison
p <- ggplot(data = q$tw) +
  geom_rect(mapping = aes(xmin = hBeg / 1000000, xmax = hEnd/1000000, 
    ymin = 0, ymax = 1, fill = hId)) +  
  layer(data = tg, geom = 'rect', mapping = 
    aes(xmin = hBeg / 1000000, xmax = hEnd / 1000000, ymin = -1, ymax = 0), 
        geom_params = list()) +
  scale_x_continuous(name = 'Chr Position (Mbp)', expand = c(0.01, 0)) +
  scale_y_continuous(name  ='', expand = c(0.04, 0)) +
  facet_grid(hId ~ .) + 
  theme(legend.position = 'right', legend.title = element_blank()) +
  theme(axis.text.x = element_text(siz e= 8, angle = 0)) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank())
ggsave(p, filename = file.path(q$dir, "figs/03_coverage.png"), 
    width = 7, height = 5)

# compare SNPs called by 2 approaches
qname = "hm340"
cq = read_comp_stat(qname, tname)
vq = read_var_stat(qname)

#fc = file.path(cq$dir, '29.gax')
#tc = read.table(fc, header = F, sep = "\t", as.is = T)[, 1:3]
#colnames(tc) = c('chr', 'beg', 'end')
#grc = GRanges(seqnames = tc$chr, ranges = IRanges(tc$beg, end = tc$end))

fs1 = file.path(cq$dir, '29.vnt/snp')
ts1 = read.table(fs1, header = F, sep = "\t", as.is = T)[, 1:2]
colnames(ts1) = c('chr', 'pos')
gs1 = GRanges(seqnames = ts1$chr, ranges = IRanges(ts1$pos, end = ts1$pos))

fs2 = file.path(vq$dir, 'snp')
ts2 = read.table(fs2, header = F, sep = "\t", as.is = T)[, 1:2]
colnames(ts2) = c('chr', 'pos')
gs2 = GRanges(seqnames = ts2$chr, ranges = IRanges(ts2$pos, end = ts2$pos))
#gs2s = intersect(grc, gr2)

sum(width(gs1))
sum(width(gs2))
sum(width(intersect(gs1, gs2)))

#fi1 = file.path(cq$dir, '29.vnt/idm')
#ti1 = read.table(fi1, header = F, sep = "\t", as.is = T)[, 1:3]
#colnames(ti1) = c('chr', 'beg', 'end')
#ti1$end = ti1$end - 1
#gi1 = GRanges(seqnames = ti1s$chr, ranges = IRanges(ti1s$beg, end = ti1s$end))

#fi2 = file.path(vq$dir, 'idm')
#ti2 = read.table(fi2, header = F, sep = "\t", as.is = T)[, 1:3]
#colnames(ti2) = c('chr', 'beg', 'end')
#gi2 = GRanges(seqnames = ti2s$chr, ranges = IRanges(ti2s$beg, end = ti2s$end))


# IDM statistics indel/other
tname = "hm101"
qname = "hm340.apecca"

dir = sprintf("/home/youngn/zhoup/Data/misc3/%s_%s/23_blat/31.9", 
  toupper(qname), toupper(tname))
fs = file.path(dir, 'idm.large.cat.tbl')
ts = read.table(fs, header = F, sep = "\t", quote = "", as.is = T)
colnames(ts) = c('tid', 'tbeg', 'tend', 'id', 'qid', 'qbeg', 'qend', 'ttype', 
  'qtype')
ts = cbind(ts, 'tlen' = ts$tend - ts$tbeg - 1, 
  'qlen' = ts$qend - ts$qbeg - 1)
tu = ts[ts$ttype == 'abs' & ts$qtype == 'abs',]
tb = ts[ts$ttype == 'pre' | ts$qtype == 'pre',]

c('unb_n' = nrow(tu), 'unb_bp_tgt' = sum(tu$tlen), 'unb_bp_qry' = sum(tu$qlen),
  'bal_n' = nrow(tb), 'bal_bp_tgt' = sum(tb$tlen), 'bal_bp_qry' = sum(tb$qlen))


# plot compstat and save to comp.png
dirstat = '/home/youngn/zhoup/Data/misc3/compstat'
tt = read.table(file.path(dirstat, 'comp.tbl'), header = T, sep = "\t", as.is = T)
tp = cbind(tt[,c(1,5)], ovl = tt$ovl / tt$ass * 100, 
  ass = (tt$ass - tt$ovl) / tt$ass * 100,
  map = (tt$map - tt$ovl) / tt$ass * 100)
tl = reshape(tp, idvar = c("sam", "type"), varying=list(3:5), timevar="prop", 
  v.names = 'pct', times = colnames(tp)[3:5], direction = 'long')

p <- ggplot(data = tl) +
  geom_bar(mapping = aes(x = sam, y = pct, fill = prop), stat = 'identity', 
    position = 'stack') +
#  layer(data=s3, geom="text", mapping=aes(x=acc, y=total/1000000 - 5, label=sprintf("%.01fx", cov_mapping)), geom_params=list(size=3.5, color='white')) +
#  layer(data=s3, geom="text", mapping=aes(x=acc, y=total/1000000 + 5, label=sprintf("%.01fx", cov_sequencing)), geom_params=list(size=3.5)) +
  scale_fill_brewer(palette = 'Set1', breaks = c("ovl", 'ass', 'map'), 
    labels=c("SNPs called by both approaches", "SNPs only called by assembly-based approach", "SNPs only called by mapping-based approach")) +
  scale_x_discrete(name = '') +
  scale_y_continuous(name = 'Proportion') +
  coord_flip() +
  theme(legend.title = element_blank(), legend.direction = 'vertical', legend.position = 'top')
ggsave(p, filename=file.path(dirstat, "compstat.png"), width=8, height=3.5)