require(rtracklayer)
require(plyr)
require(seqinr)
require(GenomicRanges))

diro = file.path(Sys.getenv("misc3"), "compstat")

tname = "HM101"
qnames = c(
  "hm004", "hm010", "hm018", "hm022", "hm034", 
  "hm050", "hm056", "hm058", "hm060", "hm095", 
  "hm125", "hm129", "hm185", "hm324", "hm340")
qnames = c("hm056", "hm056.ap", "hm340", "hm340.ap")
qnames = toupper(qnames)

###### basic assembly stats
library(Biostrings)
source(paste("http://faculty.ucr.edu/~tgirke/",
    "Documents/R_BioCond/My_R_Scripts/contigStats.R", sep=''))

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
  
  scf = file.path(dir, "11_genome.fas")
  assembly <- readDNAStringSet(scf, "fasta")
  N <- list(acc = width(assembly))
  reflength <- sapply(N, sum)
  stat <- contigStats(N = N, reflength = reflength, style = "data")
  st = as.integer(stat$Contig_Stats)
  stats[[qname]] = matrix(c(total_len, total_bases, st[c(8,2,6,4)]), 
    nrow = 1, dimnames = list(NULL, c("Total Span", "Total Bases",
    "# Scaffolds" ,"Scaffold N50" , "Scaffold Median" , "Scaffold Max")))
}
ds = do.call(rbind.data.frame, stats)
do = ds
for (i in 1:ncol(do)) { do[,i] = format(do[,i], big.mark = ",") }
fo = file.path(diro, "01_assembly_stat.tbl")
write.table(do, fo, sep = "\t", row.names = T, col.names = T, quote = F)

#### plot scaffold size distribution
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
  
  dtn = table(tg$cat2)
  x = c()
  x[fams] = 0
  x[names(dtn)] = dtn
  
  tf = file.path(dir, "51.fas")
  y = read.fasta(tf, as.string = TRUE, seqtype = "AA")
  dl = ldply(y, nchar)
  mean_prot = mean(dl$V1)
  med_prot = median(dl$V1)
  
  stats[[qname]] = matrix(c(ngene, med_prot, x),
    nrow = 1, dimnames = list(NULL, c("Total Genes", "Median Prot Length",
    fams)))
}
ds = do.call(rbind.data.frame, stats)
do = ds
for (i in 1:ncol(do)) { do[,i] = format(do[,i], big.mark = ",") }
fo = file.path(diro, "03_annotation.tbl")
write.table(do, fo, sep = "\t", row.names = T, col.names = T, quote = F)

##### NBS-LRR stats
stats = list()

fi = file.path(Sys.getenv("misc2"), "genefam", "nbs.info")
fams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
for (qname in c(tname, qnames)) {
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  fg = file.path(dir, "42.nbs", "12.gtb")
  tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,15:17)]
  dtn = table(tg$cat3)
  x = c()
  x[fams] = 0
  x[names(dtn)] = dtn
  stats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, fams))
}
ds = do.call(rbind.data.frame, stats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
for (i in 1:ncol(do)) { do[,i] = format(do[,i], big.mark = ",") }
fo = file.path(diro, "04_nbs.tbl")
write.table(do, fo, sep = "\t", row.names = T, col.names = T, quote = F)

##### CRP stats
stats = list()

fi = file.path(Sys.getenv("misc2"), "genefam", "crp.info")
fams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
for (qname in c(tname, qnames)) {
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  fg = file.path(dir, "51.gtb")
  tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,15:17)]
  tn = tg[tg$cat2 == 'crp',]
  dtn = table(tn$cat3)
  x = c()
  x[fams] = 0
  x[names(dtn)] = dtn
  stats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, fams))
}
ds = do.call(rbind.data.frame, stats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
for (i in 1:ncol(do)) { do[,i] = format(do[,i], big.mark = ",") }
fo = file.path(diro, "05_crp.tbl")
write.table(do, fo, sep = "\t", row.names = T, col.names = T, quote = F)


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

  dir = sprintf("%s/%s_%s/23_blat", Sys.getenv("misc3"), 
    toupper(qname), toupper(tname))
  fi = file.path(dir, '41.5/gax')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:3,10)]
  colnames(ti) = c('tid', 'tbeg', 'tend', 'lev')
  aligned = sum(ti$tend - ti$tbeg + 1)
  pct_aligned = aligned / total_bases * 100

  fi = file.path(dir, '31.9/gax')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:3,10)]
  colnames(ti) = c('tid', 'tbeg', 'tend', 'lev')
  synteny = sum(ti$tend - ti$tbeg + 1)
  pct_synteny = synteny / total_bases * 100

  total_bases = format(total_bases, big.mark = ",")
  brep = format(brep, big.mark = ",")
  aligned = format(aligned, big.mark = ",")
  synteny = format(synteny, big.mark = ",")
  pct_rep = sprintf("%.01f%%", pct_rep)
  pct_aligned = sprintf("%.01f%%", pct_aligned)
  pct_synteny = sprintf("%.01f%%", pct_synteny)
  stats[[qname]] = matrix(c(total_bases, brep, pct_rep, 
    aligned, pct_aligned, synteny, pct_synteny), nrow = 1, 
    dimnames = list(NULL, c("Total Bases", 
    "Repeats", "", "Alignable to HM101", "", "Synteny Blocks", "")))
}
ds = do.call(rbind.data.frame, stats)
fo = file.path(diro, "11_comp_stat.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)


##### variation stats
stats = list()
for (qname in qnames) {
  dir = sprintf("%s/%s_%s/23_blat", Sys.getenv("misc3"), 
    toupper(qname), toupper(tname))
  fi = file.path(dir, '31.9/gax')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:3,10)]
  colnames(ti) = c('tid', 'tbeg', 'tend', 'lev')
  ti = ti[ti$lev <= 2,]
  ali = sum(ti$tend - ti$tbeg + 1)
  
  fi = file.path(dir, '31.9/snp')
  ti = read.table(fi, header = F, sep = "\t", as.is = T)[,c(1:2,8)]
  colnames(ti) = c('chr', 'beg', 'lev')
  snp = sum(ti$lev <= 2)
  snpd = sprintf("%.02f%%", snp / ali * 100)

  fi = file.path(dir, '31.9/sv.ins.tbl')
  ti = read.table(fi, header = T, sep = "\t", as.is = T)[,c(1,4)]
  l1 = strsplit(ti$loc, ":")
  str2 = matrix(unlist(l1), ncol = 2, byrow = T)[,2]
  l2 = strsplit(str2, "-")
  str3 = matrix(unlist(l2), ncol = 2, byrow = T)
  ti = cbind(ti, beg = as.numeric(str3[,1]), end = as.numeric(str3[,2]))
  ti = cbind(ti, len = ti$end - ti$beg + 1)
  tis = ti[ti$len < 50,]
  si_n = nrow(tis)
  si_b = sum(tis$len)
  til = ti[ti$len >= 50,]
  li_n = nrow(til)
  li_b = sum(til$len)


  fi = file.path(dir, '31.9/sv.gan.tbl')
  ti = read.table(fi, header = T, sep = "\t", as.is = T)[,c(1,4)]
  l1 = strsplit(ti$loc, ":")
  str2 = matrix(unlist(l1), ncol = 2, byrow = T)[,2]
  l2 = strsplit(str2, "-")
  str3 = matrix(unlist(l2), ncol = 2, byrow = T)
  ti = cbind(ti, beg = as.numeric(str3[,1]), end = as.numeric(str3[,2]))
  ti = cbind(ti, len = ti$end - ti$beg + 1)
  g_n = nrow(ti)
  g_b = sum(ti$len)

  fi = file.path(dir, '31.9/sv.del.tbl')
  ti = read.table(fi, header = T, sep = "\t", as.is = T)[,c(1:3)]
  ti = cbind(ti, len = ti$end - ti$beg + 1)
  tis = ti[ti$len < 50,]
  sd_n = nrow(tis)
  sd_b = sum(tis$len)
  til = ti[ti$len >= 50,]
  ld_n = nrow(til)
  ld_b = sum(til$len)

  fi = file.path(dir, '31.9/sv.los.tbl')
  ti = read.table(fi, header = T, sep = "\t", as.is = T)[,c(1:3)]
  ti = cbind(ti, len = ti$end - ti$beg + 1)
  l_n = nrow(ti)
  l_b = sum(ti$len)
  
  fi = file.path(dir, '31.9/sv.tlc.tbl')
  ti = read.table(fi, header = T, sep = "\t", as.is = T)[,c(1:3)]
  ti = cbind(ti, len = ti$tend - ti$tbeg + 1)
  t_n = nrow(ti)
  t_b = sum(ti$len)
  
  snpstr = sprintf("%s | %5s", format(snp, big.mark = ","), snpd)
  sistr = sprintf("%s | %11s", format(si_n, big.mark = ","), format(si_b, big.mark = ","))
  sdstr = sprintf("%s | %11s", format(sd_n, big.mark = ","), format(sd_b, big.mark = ","))
  listr = sprintf("%s | %11s", format(li_n, big.mark = ","), format(li_b, big.mark = ","))
  ldstr = sprintf("%s | %11s", format(ld_n, big.mark = ","), format(ld_b, big.mark = ","))
  gstr = sprintf("%s | %11s", format(g_n, big.mark = ","), format(g_b, big.mark = ","))
  lstr = sprintf("%s | %11s", format(l_n, big.mark = ","), format(l_b, big.mark = ","))
  tstr = sprintf("%s | %11s", format(t_n, big.mark = ","), format(t_b, big.mark = ","))
  stat = c(snpstr, sdstr, sistr, ldstr, listr, lstr, gstr, tstr)
  stats[[qname]] = matrix(stat, nrow = 1, dimnames = list(NULL, 
    c("SNP (Number | Density)", "Small Del", "Small Ins", "Large Del", 
    "Large Ins", "Copy Number Loss", "Copy Number Gain", "Translocation")))
}
ds = do.call(rbind.data.frame, stats)
fo = file.path(diro, "15_vnt_stat.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)


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
  
  dir = sprintf("%s/%s_%s/41_novseq", Sys.getenv("misc3"), 
    toupper(qname), toupper(tname))
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