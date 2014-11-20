require(rtracklayer)
require(plyr)
require(seqinr)
require(GenomicRanges)

diro = file.path(Sys.getenv("misc3"), "compstat")

tname = "HM101"
qnames = c(
  "HM058", "HM125", "HM056", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340"
)
#qnames = c("HM056", "HM056.AP", "HM340", "HM340.AP")

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
  
  if(qname == tname) {
    scf_stat = rep('', 4)
    ctg_stat = rep('', 4)
  } else {
    f_scf = file.path(dir, "11_genome.fas")
    assembly <- readDNAStringSet(f_scf, "fasta")
    N <- list(acc = width(assembly))
    reflength <- sapply(N, sum)
    st <- contigStats(N = N, reflength = reflength, style = "data")
    scf_stat = as.integer(st$Contig_Stats)[c(8,2,6,4)]
    
    f_ctg = file.path(dir, "ctg.fas")
    assembly <- readDNAStringSet(f_ctg, "fasta")
    N <- list(acc = width(assembly))
    reflength <- sapply(N, sum)
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
    "Number" ,"N50", "Median" , "Max", "Number" ,"N50", "Median" , "Max", 
    "Bases", "Percent")))
}
ds = do.call(rbind.data.frame, stats)
fo = file.path(diro, "01_assembly_stat.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)

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

  dir = sprintf("%s/%s_%s/23_blat", Sys.getenv("misc3"), 
    toupper(qname), toupper(tname))
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
  dir = sprintf("%s/%s_%s", Sys.getenv("misc3"), toupper(qname), toupper(tname))
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

