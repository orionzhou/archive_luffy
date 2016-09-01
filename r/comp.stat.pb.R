require(rtracklayer)
require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
require(gtable)
source('Location.R')
source('comp.fun.R')

diro = file.path(Sys.getenv("misc2"), "pb.stat")

qnames = c("HM340", "HM340.PB", "HM340.PBBN", "HM340.PBDT", "HM340.PBBNDT", "HM340.PBDTBN")


tlen = read.table(tcfg$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(tcfg$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

##### functional annotation stat
stats = list()
statd = list()
nstats = list()
rstats = list()
cstats = list()
tstats = list()

orgs_rnaseq = c("HM034", "HM056", "HM101", "HM340")

fi = file.path(Sys.getenv("misc2"), "genefam", "nbs.info")
nfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("misc2"), "genefam", "crp.info")
cfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("misc2"), "genefam", "rlk.info")
rfams = read.table(fi, header = T, sep = "\t", as.is = T)[,1]
fi = file.path(Sys.getenv("data"), 'db', 'pfam', 'genefam.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
tfams = ti$dom[ti$fam == 'TE']
fams = unique(ti$fam[order(ti$pri)])

for (qname in c(tname, qnames_15)) {
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  fg = file.path(dir, "51.gtb")
  tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,15:17)]
  ngene = nrow(tg)
  ids_nonte = tg$id[tg$cat2 != "TE"]
  n_nonte = length(ids_nonte)
  
  dtn = table(tg$cat2)
  x = c()
  x[fams] = 0
  x[names(dtn)] = dtn
  statd[[qname]] = matrix(c(x, ngene),
    nrow = 1, dimnames = list(NULL, c(fams, "Total Genes")))
  
  tf = file.path(dir, "51.fas")
  y = read.fasta(tf, as.string = TRUE, seqtype = "AA")
  dl = ldply(y, nchar)
  lens = dl$V1[dl$.id %in% ids_nonte]
  mean_prot = mean(lens)
  med_prot = median(lens)

  tt = tg[tg$cat2 == 'TE',]
  dtn = table(tt$cat3)
  dtn = dtn[names(dtn) != '']
  x = c()
  x[tfams] = 0
  x[names(dtn)] = dtn
  tstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, tfams))
  n_te = nrow(tt)
  
  tn = tg[tg$cat2 %in% c('CC-NBS-LRR', 'TIR-NBS-LRR'),]
  dtn = table(tn$cat3)
  x = c()
  x[nfams] = 0
  x[names(dtn)] = dtn
  nstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, nfams))
  n_nbs = nrow(tn)
  
  tr = tg[tg$cat2 %in% c('RLK','LRR-RLK'),]
  dtn = table(tr$cat3)
  x = c()
  x[rfams] = 0
  x[names(dtn)] = dtn
  rstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, rfams))
  n_rlk = nrow(tr)
  
  tc = tg[tg$cat2 %in% c('CRP0000-1030','NCR','CRP1600-6250'),]
  dtn = table(tc$cat3)
  x = c()
  x[cfams] = 0
  x[names(dtn)] = dtn
  cstats[[qname]] = matrix(x, nrow = 1, dimnames = list(NULL, cfams))
  n_ncr = sum(tc$cat2 == 'NCR')
  
  tf = tg[tg$cat2 %in% c('F-box'),]
  n_fbx = nrow(tf)
  
  ids = ids_nonte
  ids_rna = c(); ids_hom = c()
  if(qname %in% orgs_rnaseq) {
    org = qname
    ff = file.path(Sys.getenv("misc2"), "rnaseq/mt/31_cufflinks", org, 'isoforms.fpkm_tracking')
    tf = read.table(ff, header = T, sep = "\t", as.is = T)
    ids_rna = tf$tracking_id[tf$tracking_id %in% ids & tf$FPKM > 0]
    p_rna = sprintf("%.01f", length(ids_rna) / length(ids) * 100)
  } else {
    p_rna = '-'
  }
  if(qname != tname) {
    ft = sprintf("%s/%s_HM101/51_ortho/31.ortho.tbl", Sys.getenv("misc3"), qname)
    tt = read.table(ft, header = T, sep = "\t", as.is = T)
    ids_hom = tt$qid[tt$qid %in% ids]
    p_hom = sprintf("%.01f", length(ids_hom) / length(ids) * 100)
  } else {
    p_hom = '-'
  }
  p_sup = sprintf("%.01f", length(unique(c(ids_rna, ids_hom))) / length(ids) * 100)
  
  stats[[qname]] = matrix(c(ngene, n_te, n_nonte, n_nbs, n_fbx, n_rlk, n_ncr, med_prot, p_rna, p_hom, p_sup), nrow = 1, dimnames = list(NULL, c("# Total Genes", 'TE', 'non-TE','NBS-LRR','F-box','LRR-RLK','NCR',"Median Prot Length", 'RNA-seq (%)', 'Homology (%)', 'RNA-seq + Homology (%)')))
}
ds = do.call(rbind.data.frame, stats)
#for (i in 1:ncol(do)) { do[,i] = format(do[,i], big.mark = ",") }
fo = file.path(diro, "03_annotation.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)

ds = do.call(rbind.data.frame, statd)
fo = file.path(diro, "04.annotation.detail.tbl")
write.table(t(ds), fo, sep = "\t", row.names = T, col.names = T, quote = F)

ds = do.call(rbind.data.frame, nstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "04.nbs.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

ds = do.call(rbind.data.frame, rstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "04.rlk.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

ds = do.call(rbind.data.frame, cstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "04.crp.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

ds = do.call(rbind.data.frame, tstats)
colsums = apply(ds, 2, sum)
do = ds[,colsums>0]
do = t(do)
do = cbind('sub-family' = rownames(do), do)
fo = file.path(diro, "04.te.tbl")
write.table(do, fo, sep = "\t", row.names = F, col.names = T, quote = F)

##### comparative / novel-seq stats
stats = list()
for (qname in qnames) {
  #add repeatmasker stats
  dir = sprintf("%s/%s", Sys.getenv("genome"), qname)
  flen = file.path(dir, "15.sizes")
  tlen = read.table(flen, sep = "\t", header = F, as.is = T)
  total_len = sum(tlen$V2)
  
  fgap = file.path(dir, "16.gap.bed")
  if(file.info(fgap)$size == 0) {
  	tgap = data.frame(V1=c(),V2=c(),V3=c())
  } else {
  	tgap = read.table(fgap, sep = "\t", header = F, as.is = T)
  }
  total_gap = sum(tgap$V3 - tgap$V2)
  total_bases = total_len - total_gap
  
  fcds = file.path(dir, "51.tbl")
  #tcds = read.table(fcds, sep = "\t", header = F, as.is = T)[,1:6]
  #tcds = tcds[tcds$V6 == 'cds',]
  #grc = GRanges(seqnames = tcds$V1, ranges = IRanges(tcds$V2, end = tcds$V3))
  #grc = reduce(grc)
  
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
  ti = ti[ti$lev <= 2,]
  synteny = sum(ti$tend - ti$tbeg + 1)
#  pct_synteny = synteny / total_bases * 100

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
  brep = format(brep, big.mark = ",")
#  pct_rep = sprintf("%.01f%%", pct_rep)
  aligned = format(aligned, big.mark = ",")
#  pct_aligned = sprintf("%.01f%%", pct_aligned)
  synteny = format(synteny, big.mark = ",")
#  pct_synteny = sprintf("%.01f%%", pct_synteny)
  bnov = format(bnov, big.mark = ",")
  pnov = sprintf("%.01f%%", pnov)
  bcds = 0#format(bcds, big.mark = ",")
  pcds = 0#sprintf("%.01f%%", pcds)
  stats[[qname]] = matrix(c(total_bases, brep, aligned, synteny, bnov, pnov, bcds, pcds), 
    nrow = 1, dimnames = list(NULL, c("Total Bases", "Repetitive", "Alignable to HM101", "Bases in Synteny Blocks", "Novel Sequences", "", "Novel Coding Seq", "")))
}
ds = do.call(rbind.data.frame, stats)
fo = file.path(diro, "11_comp_stat.tbl")
write.table(ds, fo, sep = "\t", row.names = T, col.names = T, quote = F)

