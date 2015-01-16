require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(gtable)
#require(Gviz)
source('Location.R')
source('comp.fun.R')

diro = file.path(Sys.getenv("misc3"), "compstat")

tname = "HM101"
qnames = c(
  "HM058", "HM125", "HM056", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340"
)
#qnames = c("HM056", "HM056.AC", "HM340", "HM340.AC")
#qnames = c("HM056", "HM324", "Malbus")

###### basic assembly stats
library(Biostrings)
source(paste("http://faculty.ucr.edu/~tgirke/",
    "Documents/R_BioCond/My_R_Scripts/contigStats.R", sep=''))

reflength <- 500000000
stats = list()
for (qname in c(tname, qnames)) {
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  flen = file.path(dir, "15.sizes")
  tlen = read.table(flen, sep = "\t", header = F, as.is = T)
  total_len = sum(tlen$V2)
  
  fgap = file.path(dir, "16.gap.bed")
  tgap = read.table(fgap, sep = "\t", header = F, as.is = T)
  total_gap = sum(tgap$V3 - tgap$V2)
  total_bases = total_len - total_gap
  
  if(qname == tname) {
    scf_stat = rep(0, 4)
    ctg_stat = rep(0, 4)
  } else {
    f_scf = file.path(dir, "11_genome.fas")
    assembly <- readDNAStringSet(f_scf, "fasta")
    N <- list(acc = width(assembly))
#    reflength <- sapply(N, sum)
    st <- contigStats(N = N, reflength = reflength, style = "data")
    scf_stat = as.integer(st$Contig_Stats)[c(8,2,6,4)]
    
    f_ctg = file.path(dir, "ctg.fas")
    assembly <- readDNAStringSet(f_ctg, "fasta")
    N <- list(acc = width(assembly))
#    reflength <- sapply(N, sum)
    st <- contigStats(N = N, reflength = reflength, style = "data")
    ctg_stat = as.integer(st$Contig_Stats)[c(8,2,6,4)]
  }
  
  # repeat stats
  frep = file.path(dir, "12.rm.bed")
  trep = read.table(frep, sep = "\t", header = F, as.is = T)
  brep = sum(trep$V3 - trep$V2)
  pct_rep = brep / total_bases * 100
  pct_rep = sprintf("%.02f", pct_rep)

  stat = c(total_len, total_bases, scf_stat, ctg_stat, brep)
  stat = format(stat, big.mark = ",")
  stat = c(stat, pct_rep)
  stats[[qname]] = matrix(stat, nrow = 1, 
    dimnames = list(NULL, c("Total Span", "Total Bases",
    "Number", "N50", "Median", "Max", "Number", "N50", "Median", "Max", 
    "Bases", "Percent")))
}
ds = do.call(rbind.data.frame, stats)
fo = file.path(diro, "01_assembly_stat.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)

### plot scaffold size distribution
tmp = cut(t_len$length / 1000, breaks = c(0,1,5,10,50,100,500,1000,5000))
p = ggplot(data.frame(size = tmp)) +
  geom_bar(aes(x = factor(size)), width = 0.7) + 
  scale_x_discrete(name = "Scaffold Size (kb)") + 
  scale_y_continuous(name = "") +
  theme(axis.text.x = element_text(angle = 15, size = 8))
ggsave(file.path(q$dir, "figs/01_scaffold_size.png"), p, width=5, height=4)


##### functional annotation stats
stats = list()

fi = file.path(Sys.getenv("data"), "db", "pfam", 'genefam.tbl')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
fams = unique(ti$fam[order(ti$pri)])
for (qname in c(tname, qnames)) {
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  fg = file.path(dir, "51.gtb")
  tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,15:17)]
  ngene = nrow(tg)
  ids_nonte = tg$id[tg$cat2 != "TE"]
  
  dtn = table(tg$cat2)
  x = c()
  x[fams] = 0
  x[names(dtn)] = dtn
  
  tf = file.path(dir, "51.fas")
  y = read.fasta(tf, as.string = TRUE, seqtype = "AA")
  dl = ldply(y, nchar)
  lens = dl$V1[dl$.id %in% ids_nonte]
  mean_prot = mean(lens)
  med_prot = median(lens)
  
  stats[[qname]] = matrix(c(ngene, med_prot, x),
    nrow = 1, dimnames = list(NULL, c("Total Genes", "Median Prot Length",
    fams)))
}
ds = do.call(rbind.data.frame, stats)
do = t(ds)
#for (i in 1:ncol(do)) { do[,i] = format(do[,i], big.mark = ",") }
fo = file.path(diro, "03_annotation.tbl")
write.table(do, fo, sep = "\t", row.names = T, col.names = T, quote = F)

##### NBS-LRR / CRP / TE stats
nstats = list()
cstats = list()
tstats = list()

fi = file.path(Sys.getenv("misc2"), "genefam", "nbs.info")
nfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("misc2"), "genefam", "crp.info")
cfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("data"), 'db', 'pfam', 'genefam.tbl')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
tfams = ti$dom[ti$fam == 'TE']

for (qname in c(tname, qnames)) {
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  fg = file.path(dir, "51.gtb")
  tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,15:17)]
  
  tn = tg[tg$cat2 == 'NBS-LRR',]
  dtn = table(tn$cat3)
  x = c()
  x[nfams] = 0
  x[names(dtn)] = dtn
  nstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, nfams))
  
  tn = tg[tg$cat2 == 'CRP',]
  dtn = table(tn$cat3)
  x = c()
  x[cfams] = 0
  x[names(dtn)] = dtn
  cstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, cfams))
  
  tn = tg[tg$cat2 == 'TE',]
  dtn = table(tn$cat3)
  x = c()
  x[tfams] = 0
  x[names(dtn)] = dtn
  tstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, tfams))
}

ds = do.call(rbind.data.frame, nstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "04_nbs.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

ds = do.call(rbind.data.frame, cstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "05_crp.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

ds = do.call(rbind.data.frame, tstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "06_te.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

##### comp stats
stats = list()
for (qname in qnames) {
  
  #add repeatmasker stats
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  flen = file.path(dir, "15.sizes")
  tlen = read.table(flen, sep = "\t", header = F, as.is = T)
  total_len = sum(tlen$V2)
  
  fgap = file.path(dir, "16.gap.bed")
  tgap = read.table(fgap, sep = "\t", header = F, as.is = T)
  total_gap = sum(tgap$V3 - tgap$V2)
  total_bases = total_len - total_gap
  
  frep = file.path(dir, "12.rm.bed")
  trep = read.table(frep, sep = "\t", header = F, as.is = T)
  brep = sum(trep$V3 - trep$V2)
  pct_rep = brep / total_bases * 100

  dir = sprintf("%s/%s_%s/23_blat", Sys.getenv("misc3"), qname, tname)
  fi = file.path(dir, '41.5/gax')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:3,10)]
  colnames(ti) = c('tid', 'tbeg', 'tend', 'lev')
  aligned = sum(ti$tend - ti$tbeg + 1)
#  pct_aligned = aligned / total_bases * 100

  fi = file.path(dir, '31.9/gax')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:3,10)]
  colnames(ti) = c('tid', 'tbeg', 'tend', 'lev')
  synteny = sum(ti$tend - ti$tbeg + 1)
#  pct_synteny = synteny / total_bases * 100

  total_bases = format(total_bases, big.mark = ",")
  brep = format(brep, big.mark = ",")
  aligned = format(aligned, big.mark = ",")
  synteny = format(synteny, big.mark = ",")
  pct_rep = sprintf("%.01f%%", pct_rep)
#  pct_aligned = sprintf("%.01f%%", pct_aligned)
#  pct_synteny = sprintf("%.01f%%", pct_synteny)
  stats[[qname]] = matrix(c(total_bases, brep, pct_rep, aligned, synteny), 
    nrow = 1, dimnames = list(NULL, c("Total Bases", 
    "Repeats", "", "Alignable to HM101", "Synteny Blocks")))
}
ds = do.call(rbind.data.frame, stats)
fo = file.path(diro, "11_comp_stat.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)


##### variation stats
stats = list()
for (qname in qnames) {
  dir = sprintf("%s/%s_%s", Sys.getenv("misc3"), qname, tname)
  fi = file.path(dir, '23_blat/31.9/gax')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:3,10)]
  colnames(ti) = c('tid', 'tbeg', 'tend', 'lev')
  ti = ti[ti$lev <= 2,]
  ali = sum(ti$tend - ti$tbeg + 1)
  
  fi = file.path(dir, '23_blat/31.9/snp')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:2,8)]
  colnames(ti) = c('chr', 'beg', 'lev')
  snp = sum(ti$lev <= 2)
  snpd = sprintf("%.02f%%", snp / ali * 100)

  fv = file.path(dir, '31_sv/02.sort.stb')
  tv = read.table(fv, header = T, sep = "\t", as.is = T)

  ti = tv[tv$type == 'INS',]
  ti = cbind(ti, len = ti$qend - ti$qbeg + 1)
  tis = ti[ti$len < 50,]
  si_n = nrow(tis)
  si_b = sum(tis$len)
  til = ti[ti$len >= 50,]
  li_n = nrow(til)
  li_b = sum(til$len)

  ti = tv[tv$type == 'CNG',]
  ti = cbind(ti, len = ti$qend - ti$qbeg + 1)
  g_n = nrow(ti)
  g_b = sum(ti$len)

  ti = tv[tv$type == 'DEL',]
  ti = cbind(ti, len = ti$tend - ti$tbeg + 1)
  tis = ti[ti$len < 50,]
  sd_n = nrow(tis)
  sd_b = sum(tis$len)
  til = ti[ti$len >= 50,]
  ld_n = nrow(til)
  ld_b = sum(til$len)

  ti = tv[tv$type == 'CNL',]
  ti = cbind(ti, len = ti$tend - ti$tbeg + 1)
  l_n = nrow(ti)
  l_b = sum(ti$len)
  
  fi = file.path(dir, '31_sv/05.refine.tlc')
  ti = read.table(fi, header = T, sep = "\t", as.is = T)[,c(1:3)]
  ti = cbind(ti, len = ti$tdend - ti$tdbeg + 1)
  t_n = nrow(ti)
  t_b = sum(ti$len)
  
  stat = c(snp, si_n, si_b, sd_n, sd_b, li_n, li_b, ld_n, ld_b,
    g_n, g_b, l_n, l_b, t_n, t_b)
  stat = format(stat, big.mark = ",")
  stat = c(stat[1], snpd, stat[-1])
  stats[[qname]] = matrix(stat, nrow = 1, dimnames = list(NULL, 
    c("SNP #", "SNP Density (/bp)", "Small Del (#)", "Small Del (bp)", 
    "Small Ins (#)", "Small Ins (bp)", "Large Del (#)", "Large Del (bp)", 
    "Large Ins (#)", "Large Ins (bp)", 
    "Copy Number Loss (#)", "Copy Number Loss (bp)", 
    "Copy Number Gain (#)", "Copy Number Gain (bp)", 
    "Translocation (#)", "Translocation (bp)")))
}
ds = do.call(rbind.data.frame, stats)
fo = file.path(diro, "15_vnt_stat.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)


##### create sliding window table
tcfg = get_genome_cfg(tname)
tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

x = tt$end
names(x) = tt$chr
gr = tileGenome(x, tilewidth = 100000, cut.last.tile.in.chrom = T)

tw = data.frame(chr = seqnames(gr), beg = start(gr), end = end(gr), 
  len = width(gr), stringsAsFactors = F)
  
bp_gap = intersect_basepair(gr, grp)
bp_nogap = tw$len - bp_gap

to = cbind(tw, len_ng = bp_nogap)
fo = file.path(diro, "31.win.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### test correlation of pi ~ gene density
fw = file.path(Sys.getenv("misc3"), "compstat", "32.win.stat.tbl")
tw = read.table(fw, header = T, sep = "\t", as.is = T)
gr = with(tw, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

tg = read.table(tcfg$gene, sep = "\t", header = F, as.is = T)[,1:7]
colnames(tg) = c("chr", 'beg', 'end', 'srd', 'id', 'type', 'cat')
tg = tg[tg$type == 'cds',]
tgg = tg[tg$cat != 'TE',]
tgt = tg[tg$cat == 'TE',]
tgn = tg[tg$cat == 'NBS-LRR',]
tgc = tg[tg$cat == 'CRP',]
ggg = with(tgg, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
ggt = with(tgt, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
ggn = with(tgn, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
ggc = with(tgc, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

p_g = intersect_basepair(gr, reduce(ggg)) / tw$len_ng
p_t = intersect_basepair(gr, reduce(ggt)) / tw$len_ng
p_n = intersect_basepair(gr, reduce(ggn)) / tw$len_ng
p_c = intersect_basepair(gr, reduce(ggc)) / tw$len_ng
tgd = data.frame('gene' = p_g, 'te' = p_t, 'nbs' = p_n, 'crp' = p_c, as.is = T)

fit = lm(tw$pi ~ p_g + p_t + p_n + p_c)
summary(fit)

# test correlation of gene density with distance to centromere
ts = tw[tw$chr == 'chr5',]
to = transform(ts, pct_cds = bp_cds/bp_nogap, bp_dist = abs( (beg+end)/2 - 21450000))
fit = lm(pct_cds ~ bp_dist, data = to)
summary(fit)
plot(to$bp_dist, to$pct_cds)

# test non-NCR v.s. NCR
fb = file.path(Sys.getenv("genome"), tname, "51.gtb")
tb = read.table(fb, sep = "\t", header = T, as.is = T)[,c(1,17)]
tg2 = merge(tg, tb, by = 'id')

tgc1 = tg2[tg2$cat == 'CRP' & tg2$cat3 <= 'CRP1030',]
tgc2 = tg2[tg2$cat == 'CRP' & tg2$cat3 > 'CRP1030' & tg2$cat3 < 'CRP1600',]
tgc3 = tg2[tg2$cat == 'CRP' & tg2$cat3 >= 'CRP1600',]
ggc1 = with(tgc1, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
ggc2 = with(tgc2, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
ggc3 = with(tgc3, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

p_c1 = intersect_basepair(gr, reduce(ggc1)) / bp_nogap
p_c2 = intersect_basepair(gr, reduce(ggc2)) / bp_nogap
p_c3 = intersect_basepair(gr, reduce(ggc3)) / bp_nogap

pc = p_c1 + p_c2 + p_c3
fit <- lm(tw$pi ~ tw$gen + tw$tre + tw$nbs + p_c2)
summary(fit)

##### variant density heatmap
fw = file.path(Sys.getenv("misc3"), "compstat", "32.win.stat.tbl")
tw = read.table(fw, header = T, sep = "\t", as.is = T)
gr = with(tw, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

sdel = list(); ldel = list(); sins = list(); lins = list()
cng = list(); cnl = list(); tlc = list()
qnames = c("HM058", "HM034", "HM022")
for (qname in qnames) {
  dir = sprintf("%s/%s_%s", Sys.getenv("misc3"), qname, tname)
  
  fv = file.path(dir, '31_sv/02.sort.stb')
  tv = read.table(fv, header = T, sep = "\t", as.is = T)

  ti = tv[tv$type == 'DEL',]
  ti = cbind(ti, len = ti$tend - ti$tbeg + 1)
  tis = ti[ti$len < 50,]
  til = ti[ti$len >= 50,]
  grs = with(tis, GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend)))
  grl = with(til, GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend)))
  bp_sdel = intersect_basepair(gr, grs)
  bp_ldel = intersect_basepair(gr, grl)
  sdel[[qname]] = bp_sdel / tw$len_ng
  ldel[[qname]] = bp_ldel / tw$len_ng
  
  ti = tv[tv$type == 'INS',]
  ti = cbind(ti, len = ti$qend - ti$qbeg + 1)
  tis = ti[ti$len < 50,]
  til = ti[ti$len >= 50,]
  grs = with(tis, GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend), score = len))
  grl = with(til, GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend), score = len))
  bp_sins = intersect_score(gr, grs)
  bp_lins = intersect_score(gr, grl)
  sins[[qname]] = bp_sins / tw$len_ng
  lins[[qname]] = bp_lins / tw$len_ng
  
  ti = tv[tv$type == 'CNL',]
  ti = cbind(ti, len = ti$tend - ti$tbeg + 1)
  grs = with(ti, GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend)))
  bp_cnl = intersect_basepair(gr, grs)
  cnl[[qname]] = bp_cnl / tw$len_ng

  ti = tv[tv$type == 'CNG',]
  ti = cbind(ti, len = ti$qend - ti$qbeg + 1)
  grs = with(ti, GRanges(seqnames = tchr, ranges = IRanges(tbeg, end = tend), score = len))
  bp_cng = intersect_score(gr, grs)
  cng[[qname]] = bp_cng / tw$len_ng
  
  fv = file.path(dir, '31_sv/05.refine.tlc')
  ti = read.table(fv, header = T, sep = "\t", as.is = T)
  grs = with(ti, GRanges(seqnames = tdchr, ranges = IRanges(tdbeg, end = tdend)))
  bp_tlc = intersect_basepair(gr, grs)
  tlc[[qname]] = bp_tlc / tw$len_ng
}
sdel = do.call(rbind, sdel); ldel = do.call(rbind, ldel)
sins = do.call(rbind, sins); lins = do.call(rbind, lins)
cng = do.call(rbind, cng); cnl = do.call(rbind, cnl); tlc = do.call(rbind, tlc)

sid = sins + sdel; lid = lins + ldel; cnv = cng + cnl


chr = 'chr5'
idxs = which(tw$chr == chr)
to = tw[idxs,]
to$pi[to$lenc < 2000] = NA

gd = tgd[idxs,]

labs = seq(0, floor(tt$end[tt$chr == chr]) / 1000000, by = 10)
pb <- ggplot(cbind(to, gd)) +
  theme_bw() + 
  scale_x_continuous(name = 'chromosome position (Mbp)', expand = c(0, 0), breaks = labs*1000000, labels = labs) +
  theme(axis.title.x = element_text(colour = "black", size = 9)) +
  theme(axis.text.x = element_text(colour = "blue", size = 8)) +
  theme(axis.title.y = element_text(colour = "blue", size = 9)) +
  theme(axis.text.y = element_text(colour = "grey", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))

p_ng <- pb +
  geom_rect(aes(xmin = beg, xmax = end, ymin = 0, ymax = len_ng/len)) +
  scale_y_continuous(expand = c(0, 0), name = 'Gap-free bases')
p_gd1 <- pb +
  geom_rect(aes(xmin = beg, xmax = end, ymin = 0, ymax = gene, fill = 'Gene')) +
  geom_rect(aes(xmin = beg, xmax = end, ymin = gene, ymax = gene+te, fill = 'TE')) +
  scale_y_continuous(expand = c(0, 0), name = 'Gene density')
p_gd2 <- pb +
  geom_rect(aes(xmin = beg, xmax = end, ymin = 0, ymax = nbs, fill = 'NBS-LRR')) +
  geom_rect(aes(xmin = beg, xmax = end, ymin = nbs, ymax = nbs+crp, fill = 'CRP')) +
  scale_y_continuous(expand = c(0, 0), name = 'Gene density')
p_c <- pb +
  geom_rect(aes(xmin = beg, xmax = end, ymin = 0, ymax = lenc/len)) +
  scale_y_continuous(expand = c(0, 0), name = 'Covered bases')
p_pi <- pb +
  geom_rect(aes(xmin = beg, xmax = end, ymin = 0, ymax = pi)) +
  scale_y_continuous(expand = c(0, 0), name = 'Pi')

gt_ng <- ggplot_gtable(ggplot_build(p_ng)) 
gt_gd1 <- ggplot_gtable(ggplot_build(p_gd1)) 
gt_gd2 <- ggplot_gtable(ggplot_build(p_gd2)) 
gt_c <- ggplot_gtable(ggplot_build(p_c)) 
gt_pi <- ggplot_gtable(ggplot_build(p_pi)) 

tracks = list(gt_ng, gt_gd1, gt_gd2, gt_c, gt_pi)
trackheights = c(5, 5, 5, 5, 5)

pad = mean(trackheights) * 0.05
hts = as.vector(rbind(rep(pad, length(tracks)), trackheights, rep(pad, length(tracks))))
gt <- gtable(widths = unit(c(1, 10, 2), "null"), height = unit(c(hts, 2), "null")) 

for (i in 1:length(tracks)) {
  gt1 = tracks[[i]]
  gt <- gtable_add_grob(gt, gt1[3, 2:3], i*3-1, 1)
  gt <- gtable_add_grob(gt, gt1[3, 4], i*3-1, 2) 
  if(ncol(gt1) == 6) {
    gt <- gtable_add_grob(gt, gt1[3, 5], i*3-1, 3)
  }
}
gt <- gtable_add_grob(gt, gt1[4:5, 4], length(tracks)*3 + 1, 2)

fo = sprintf("%s/gd.pdf", diro)
pdf(file = fo, width = 8, height = 6, bg = 'transparent')
grid.newpage() 
grid.draw(gt)
dev.off()

## Gviz plot (obsolete)
options(ucscChromosomeNames = FALSE)
axisTrack <- GenomeAxisTrack(cex = 0.8, exponent = 6, labelPos = "below")
nogapTrack <- DataTrack(chromosome = chr, start = to$beg, end = to$end, 
  genome = "Mt",  name = "nogap", data = nogap/to$len, type = "histogram")
gdTrack1 <- DataTrack(chromosome = chr, start = to$beg, end = to$end, 
  genome = "Mt",  name = "Gene/TE", data = gd1, groups = cats1, 
  col = c('palegreen', 'royalblue2'), col.histogram = 'red', lwd = 0,
  type = "histogram", stackedBars = F, legend = F, cex.legend = 0.7)
gdTrack2 <- DataTrack(chromosome = chr, start = to$beg, end = to$end, 
  genome = "Mt",  name = "NBS/CRP", data = gd2, groups = cats2, 
  col = c('orangered', 'royalblue2'), col.histogram = 'red', lwd = 0,
  type = "histogram", stackedBars = F, legend = F, cex.legend = 0.7)

gaxTrack <- DataTrack(chromosome = chr, start = to$beg, end = to$end, 
  genome = "Mt", name = "gax", data = to$lenc, type = 'histogram',
  lwd = 0.5, legend = T, cex.legend = 0.7)
piTrack <- DataTrack(chromosome = chr, start = to$beg, end = to$end, 
  genome = "Mt", name = "Pi", data = to$pi, type = 'histogram',
  lwd = 0.5)
sidTrack <- DataTrack(chromosome = chr, start = to$beg, end = to$end, 
  genome = "Mt", name = "sid", data = sid[,idxs], groups = qnames, type = 'a',
  lwd = 0.5)
lidTrack <- DataTrack(chromosome = chr, start = to$beg, end = to$end, 
  genome = "Mt", name = "lid", data = lid[,idxs], groups = qnames, type = 'a',
  lwd = 0.5)
cnvTrack <- DataTrack(chromosome = chr, start = to$beg, end = to$end, 
  genome = "Mt", name = "cnv", data = cnv[,idxs], groups = qnames, type = 'a',
  lwd = 0.5)
tlcTrack <- DataTrack(chromosome = chr, start = to$beg, end = to$end, 
  genome = "Mt", name = "tlc", data = tlc[,idxs], groups = qnames, type = 'a',
  lwd = 0.5)

fo = sprintf("%s/gd.pdf", diro)
pdf(file = fo, width = 8, height = 6, bg = 'transparent')
plotTracks(
  list(axisTrack, nogapTrack, gdTrack1, gdTrack2,
    gaxTrack, piTrack, sidTrack, lidTrack, cnvTrack, tlcTrack),
  chromosome = chr, from = 1, to = tt$end[tt$chr == chr],
  sizes = c(0.5, 0.7, 0.8, 0.8,
    0.7, 1, 1, 1, 1, 1)
)
dev.off()

##### NovelSeq stats
stats = list()
for (qname in qnames) {
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  flen = file.path(dir, "15.sizes")
  tlen = read.table(flen, sep = "\t", header = F, as.is = T)
  total_len = sum(tlen$V2)
  
  fgap = file.path(dir, "16.gap.bed")
  tgap = read.table(fgap, sep = "\t", header = F, as.is = T)
  total_gap = sum(tgap$V3 - tgap$V2)
  total_bases = total_len - total_gap
  
  fcds = file.path(dir, "51.tbl")
  tcds = read.table(fcds, sep = "\t", header = F, as.is = T)[,1:6]
  tcds = tcds[tcds$V6 == 'cds',]
  grc = GRanges(seqnames = tcds$V1, ranges = IRanges(tcds$V2, end = tcds$V3))
  grc = reduce(grc)
  
  dir = sprintf("%s/%s_%s/41_novseq", Sys.getenv("misc3"), qname, tname)
  fi = file.path(dir, '21.bed')
  ti = read.table(fi, sep = "\t", header = F, as.is = T)
  bnov = sum(ti$V3 - ti$V2)
  pnov = bnov / total_bases * 100
  
  grn = GRanges(seqnames = ti$V1, ranges = IRanges(ti$V2+1, end = ti$V3))
  grn = reduce(grn)
  bcds = sum(width(intersect(grc, grn)))
  pcds = bcds / bnov * 100

  total_bases = format(total_bases, big.mark = ",")
  bnov = format(bnov, big.mark = ",")
  pnov = sprintf("%.01f%%", pnov)
  bcds = format(bcds, big.mark = ",")
  pcds = sprintf("%.01f%%", pcds)
  stats[[qname]] = matrix(c(total_bases, bnov, pnov, bcds, pcds
    ), nrow = 1, 
    dimnames = list(NULL, c("Total Bases", 
    "Novel Sequences", "", "Novel Coding Seq", "")))
}
ds = do.call(rbind.data.frame, stats)
fo = file.path(diro, "21_novseq.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)


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

