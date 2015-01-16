require(ape)
require(seqinr)
require(rbamtools)
require(plyr)
require(Gviz)
require(Biostrings)
source("Align.R")
source("Location.R")
source("comp.fun.R")

qnames = c(
  "HM058", "HM125", "HM056.AC", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340.AC"
)
tname = "HM101"
qname = "HM340.AC"
tcfg = get_genome_cfg(tname)
qcfg = get_genome_cfg(qname)

dir = file.path(Sys.getenv("misc3"), "comp.expand", qname)

tq = read.table(qcfg$gene, header = F, sep = "\t", as.is = T)
colnames(tq) = c('chr', 'beg', 'end', 'srd', 'id', 'type', 'cat')
tqf = tq[tq$type == 'mrna', c(1:5)]

qseqlen = read.table(qcfg$size, header = F, sep = "\t",  as.is = T, col.names = c("id", "size"))
qseqinfo = read_seqinfo(qcfg$size)

fi = file.path(Sys.getenv("misc3"), "comp.ortho", "33.ortho.cat.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti[is.na(ti)] = ''

##### create CRP physical clusters
ids = ti[ti$cat2 == 'CRP', qname]
ids = ids[ids != '']
td = tqf[tqf$id %in% ids,]

chrs = sort(unique(td$chr))
chrmap = 1:length(chrs)
names(chrmap) = chrs
gpos = chrmap[td$chr] * 100000000 + as.integer((td$beg + td$end) / 2)
clus = locCluster(gpos, 15000)
td = cbind(td, clu = clus$cluster, stringsAsFactors = F)

fd = file.path(dir, "01.cluster.tbl")
write.table(td, file = fd, sep = "\t", row.names = F, col.names = T, quote = F)

##### identify CRP expansion
tfas <- read.fasta(tcfg$protein, seqtype = "AA", as.string = T, set.attributes = F)
qfas <- read.fasta(qcfg$protein, seqtype = "AA", as.string = T, set.attributes = F)

tl = ti[ti$cat2 == "CRP" & ti[,tname] != '' & ti[,qname] == '', c(tname, qname)]
tm = ti[ti$cat2 == "CRP" & ti[,tname] == '' & ti[,qname] != '', c(tname, qname)]

tm = tqf[tqf$id %in% tm[,qname],]
tm = tm[order(tm$chr, tm$beg, tm$end), ]
tm = cbind(tm, aid = NA, abeg = NA, aend = NA, bpos = NA, epos = NA, stringsAsFactors = F)
for (i in 1:nrow(tm)) {
  chr = tm$chr[i]; beg = tm$beg[i]; end = tm$end[i]; id = tm$id[i]
  pos = (beg + end) / 2
  tc = tqf[tqf$chr == chr & tqf$id != id & (tqf$beg > end | tqf$end < beg) &
    abs((tqf$beg+tqf$end)/2 - pos) <= 15000, ]
  if(nrow(tc) == 0) next
  to = data.frame(tid = id, tseq = as.character(qfas[id]), 
    qid = tc$id, qseq = as.character(qfas[tc$id]),
    abeg = tc$beg, aend = tc$end, stringsAsFactors = F)
  to = cbind(to, apos = as.integer((to$abeg + to$aend) / 2))
  ta = apply(to, 1, aa_pw_dist)
  ta = data.frame(t(ta))
  idxs = ta$mat / ta$alen >= 0.9 & min(ta$tgap/ta$len, ta$qgap/ta$len) <= 0.2
  if(sum(idxs) == 0) next
  toc = to[idxs, ]
  tos = toc[order(abs(toc$apos - pos)),]
  tm$aid[i]  = tos$qid[1]
  tm$abeg[i] = tos$abeg[1]
  tm$aend[i] = tos$aend[1]
  tm$bpos[i] = min(end, tm$aend[i])
  tm$epos[i] = max(beg, tm$abeg[i])
}
tb = tm[!is.na(tm$aid),]
fm = file.path(dir, "11.exp.tbl")
write.table(tm, file = fb, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

## find pacbio support - run "comp.expand.pl"

### plot CRP expansion - Gviz
fc = file.path(dir, "01.cluster.tbl")
tc = read.table(fc, header = T, sep = "\t", as.is = T)
fx = file.path(dir, "12.pacbio.tbl")
tx = read.table(fx, header = T, sep = "\t", as.is = T)

tx = merge(tx[tx$n_read > 0, ], tc[, c('id','clu')], by = 'id')
ty = ddply(tx, .(clu), summarise, chr = chr[1], beg = min(beg, abeg), 
  end = max(end, aend), ids = paste(id, collapse = ' '), 
  aids = paste(aid, collapse = ' '), reads = paste(reads, collapse = ' '))

dna = readDNAStringSet(qcfg$dna)

fbam = file.path(Sys.getenv("misc3"), 'pacbio', 'HM340_HM340.AC', '15.bam')
bam = bamReader(fbam, idx = T)
seqmap = getRefData(bam)

read_bam <- function(bam, seqmap, chr, beg, end) {
  refid = seqmap$ID[seqmap$SN == chr]
  reads = bamRange(bam, c(refid, beg, end))
  td <- as.data.frame(reads)
  tr = data.frame(chr = chr, beg = td$position + 1, end = NA, 
    len = nchar(td$seq), srd = "+", id = td$name, cigar = td$cigar, 
    seq = td$seq, stringsAsFactors = F)
  tr$end = tr$beg + tr$len - 1
  tr$srd[which(td$revstrand)] = "-"
  tr
}

options(ucscChromosomeNames = FALSE)
axisTrack <- GenomeAxisTrack(cex = 1.0, exponent = 3)
sTrack <- SequenceTrack(dna)

for (i in 1:nrow(ty)) {
i = 19
chr = ty$chr[i]; beg = max(1, ty$beg[i]-1000); 
end = min(ty$end[i]+1000, qseqlen$size[qseqlen$id==chr])
ids = as.character(strsplit(ty$ids[i], split = " ")[[1]])
aids = as.character(strsplit(ty$aids[i], split = " ")[[1]])
rids = unique(as.character(strsplit(ty$reads[i], split = " ")[[1]]))
gr = GRanges(seqnames = chr, ranges = IRanges(beg, end = end))

tp = read_gap(qcfg$gapz, gr)
gapTrack <- AnnotationTrack(genome = qname, 
  chromosome = tp$id, start = tp$beg, end = tp$end, name = 'gap', 
  showId = F, fill = 'maroon', background.title = "midnightblue")

tg = read_gene(qcfg$genez, gr)
tg = tg[tg$type == 'cds',]
dfg = data.frame(chromosome = tg$chr, start = tg$beg, end = tg$end,
  width = tg$end-tg$beg+1, strand = tg$srd, feature = tg$type, gene = tg$id, 
  exon = NA, transcript = tg$id, symbol = tg$id, stringsAsFactors = F)
gr.fill = rep('gray', nrow(dfg))
gr.fill[tg$id %in% aids] = 'orange'
gr.fill[tg$id %in% ids] = 'green'
grTrack <- GeneRegionTrack(dfg, genome = qname, shape = 'smallArrow',
  name = "genes", showId = T, just.group = 'below', stackHeight = 0.75,
  fill = gr.fill, fontsize = 9, max.height = 10,
  cex.group = 0.8, cex.title = 1, background.title = 'midnightblue')

tr = read_bam(bam, seqmap, chr, beg, end)
reads.col = rep('#BABABA', nrow(tr))
reads.col[tr$id %in% rids] = 'firebrick'
alTrack <- AlignmentsTrack(genome = qname, 
  chromosome = tr$chr, range = fbam, 
  name = 'PacBio', showId = F, 
  fill.reads = reads.col, lwd.reads = 0, alpha.reads = 0.8, 
  col.mismatch = 'blue', lwd.mismatch = 0, alpha.mismatch = 0.6,
  showMismatches = T, noLetters = T, max.height = 5, 
  background.title = "midnightblue", isPaired = F)

fo = sprintf("%s/fig%03d.pdf", dir, i)
pdf(file = fo, width = 8, height = 6, bg = 'transparent')
  plotTracks(
    list(axisTrack, gapTrack, grTrack, alTrack, sTrack),
    chromosome = chr, from = beg, to = end,
    min.height = 0, coverageHeight = 0.08, minCoverageHeight = 0,
    sizes = c(1, 0.4, 1, 5, 1)
  )
  dev.off()
}

bamClose(bam)

### plot CRP expansion phylogeny
fc = file.path(dir, "01.cluster.tbl")
tc = read.table(fc, header = T, sep = "\t", as.is = T)
fx = file.path(dir, "12.pacbio.tbl")
tx = read.table(fx, header = T, sep = "\t", as.is = T)

fi = file.path(Sys.getenv("misc3"), "comp.ortho", "33.ortho.cat.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti[is.na(ti)] = ''

t1 = merge(tx[tx$n_read > 0, ], tc[, c('id','clu')], by = 'id')
t2 = t1[,c('id','aid','clu')]
t3 = merge(t2, ti[,c(qname, 'cat3')], by.x = 'id', by.y = qname)

get_fam <- function(z) { tab = table(z); names(tab)[which(tab == max(tab))[1]] }
ty = ddply(t3, .(clu), summarise, fam = get_fam(cat3), 
  ids = paste(id, collapse = ' '), aids = paste(aid, collapse = ' '))

for (i in 1:nrow(ty)) {
#i = 1
fam = ty$fam[i]
ids = as.character(strsplit(ty$ids[i], split = " ")[[1]])
aids = as.character(strsplit(ty$aids[i], split = " ")[[1]])

fi = sprintf("%s/05.aln/%s.ph", dir, fam)
tree = read.tree(fi)

labs = tree$tip.label
res = strsplit(labs, split = "|", fixed = T)
orgs = sapply(res, "[", 1)
gids = sapply(res, "[", 2)

cols = rep("lightskyblue1", length(labs))
cols[orgs == qname] = "pink1"
cols[gids %in% ids] = 'firebrick'
cols2 = rep("white", length(labs))
#cols2[gids %in% aids] = 'royalblue'

ht = length(tree$tip.label) * 0.12

fo = sprintf("%s/15.phy/%03d.pdf", dir, i)
pdf(file = fo, width = 6, height = ht, bg = 'transparent')

plot(tree, font = 1, show.node.label = F, show.tip.label = F,
  x.lim = c(0, 1.5), no.margin = T, cex = 0.7)
tiplabels(labs, adj = -0.05, bg = cols, cex = 0.5)

add.scale.bar(lcol = 'black')
dev.off()
}

##### identify inserted/gained CRPs overlapping SVs
qname = 'HM340'

dirq = sprintf("%s/%s", Sys.getenv("genome"), qname)
dirt = sprintf("%s/%s", Sys.getenv("genome"), tname)
dirc = sprintf("%s/%s_%s", Sys.getenv("misc3"), qname, tname)

ft = file.path(dirt, '51.tbl')
tt = read.table(ft, header = F, sep = "\t", as.is = T)
colnames(tt) = c('chr', 'beg', 'end', 'srd', 'id', 'type', 'cat')
ttf = tt[tt$type == 'mrna', c(1:5,7)]

fq = file.path(dirq, '51.tbl')
tq = read.table(fq, header = F, sep = "\t", as.is = T)
colnames(tq) = c('chr', 'beg', 'end', 'srd', 'id', 'type', 'cat')
tqf = tq[tq$type == 'mrna', c(1:5,7)]

fr = file.path(dirc, '51_ortho', '31.ortho.tbl')
tr1 = read.table(fr, header = T, sep = "\t", as.is = T)
tr2 = merge(tr1, ttf, by.x = 'tid', by.y = 'id')
colnames(tr2)[3:7] = c('tchr', 'tbeg', 'tend', 'tsrd', 'tcat')
tr = merge(tr2, tqf, by.x = 'qid', by.y = 'id')
colnames(tr)[8:12] = c('qchr', 'qbeg', 'qend', 'qsrd', 'qcat')
tr = tr[order(tr$tchr, tr$tbeg, tr$tend), c(2:7,1,8:12)]


dirs = sprintf("%s/%s_%s/31_sv", Sys.getenv("misc3"), qname, tname)
fsv = file.path(dirs, "11.cds.tbl")
tsv = read.table(fsv, sep = "\t", as.is = T, header = T)

tp = tsv[tsv$type == 'INS' | tsv$type == 'CNG',]
tp = cbind(tp, gpct = tp$gslen / tp$glen)
tps = tp[tp$gpct >= 0.8,]
tps2 = merge(tps, tqf, by.x = 'gid', by.y = 'id')

table(tps2$cat)
tps2[tps2$cat == 'CRP',]

x=table(tps2$cat)
write.table(x, file = file.path(dirs, 'cng.tbl'), sep = "\t", row.names = F, col.names = F, quote = F)

x = tps2[tps2$cat == 'CRP',]
#x$tbeg = x$tbeg - 20000
#x$tend = x$tend + 20000
write.table(x, file = file.path(dirs, 'crp.tbl'), sep = "\t", row.names = F, col.names = T, quote = F)

##### identify crp sub-clade membership
require(geiger)
require(igraph)
source('clustertree.R')

dir = file.path(Sys.getenv("misc2"), 'genefam')
fc = file.path(dir, "11.crp.tbl")
tc = read.table(fc, sep = "\t", header = T, as.is = T)[,1:13]

to = cbind(tc, clu = NA)
fams = unique(tc$family)
for (fam in fams) {
  tcs = tc[tc$family == fam,]
  if(nrow(tcs) < 10) next
  fi = sprintf("%s/23_aln/%s.ph", dir, tolower(fam))
  ff = sprintf("%s/26_fig/%s.png", dir, tolower(fam))
  if(!file.exists(fi)) next
  
  tree = read.tree(fi)
  cls = prosperi.cluster(tree, 0.05)$membership

  ids = tree$tip.label
  aorgs = sapply(strsplit(ids, "_"), "[", 1)
  
  for (i in 1:length(ids)) {
    to$clu[to$family == fam & to$treeid == ids[i]] = cls[i]
  }

  ucl = sort(unique(cls))
  uorg = sort(unique(aorgs))

  pal = rainbow(length(ucl))
  names(pal) = as.character(ucl)

  ht = length(tree$tip.label) * 12
  png(filename = ff, width = 800, height = ht, units = 'px')
  plot(tree, font = 1, show.node.label = F, show.tip.label = F,
    label.offset = 0.01, no.margin = T, cex = 0.81)
  newtip = paste(cls, tree$tip.label, sep = ' ')
  tiplabels(newtip, adj = -0.1, bg = pal[as.character(cls)], cex = 0.8)

  add.scale.bar(lcol = 'black')
  dev.off()
}
fo = sprintf("%s/27.cl.tbl", dir)
write.table(to, file = fo, sep = "\t", row.names = F, col.names = T, quote = F)

### identify crp tandup
dir = file.path(Sys.getenv("misc2"), 'genefam')
fc = file.path(dir, "27.cl.tbl")
tc = read.table(fc, sep = "\t", header = T, as.is = T)

orts = list()
for (org in orgs) {
  ft = sprintf("%s/%s_HM101/23_blat/ortho.tbl", Sys.getenv("misc3"), org)
  tt = read.table(ft, header = T, sep = "\t", as.is = T)[,c(1,3)]
  orts[[org]] = tt
}

dat = data.frame()
fams = unique(tc$family)
fo = sprintf("%s/28.expand.tbl", dir)
for (fam in fams) {
  tc1 = tc[tc$family == fam,]
  if(nrow(tc1) < 10) next

  for (cl in sort(unique(tc1$cl))) {
    if(cl == 0) next
    tc2 = tc1[tc1$cl == cl,]
  
    tcr = tc2[tc2$org == "HM101",]
    if(nrow(tcr) == 0) next
    for (org in sort(unique(tc2$org))) {
      if(org == 'HM101') next
      tc3 = tc2[tc2$org == org,]
      
      ort = orts[[org]]
      
      t_ort = ort[ort$qid %in% tc3$id,]
      ids_orphan = t_ort$qid[t_ort$tid == ""]
      if(length(ids_orphan) == 0) next
      tc4 = tc3[tc3$id %in% ids_orphan,]
      for (i in 1:nrow(tcr)) {
        idr = tcr$id[i]
        ido = ort$qid[ort$tid == idr]
        chr = tc3$chr[tc3$id == ido]
        pos = (tc3$beg[tc3$id == ido] + tc3$end[tc3$id == ido]) / 2
        
        if(ido == '') next
        idxs = tc4$chr == chr & abs((tc4$beg+tc4$end)/2 - pos) < 5000
        if(sum(idxs) == 0) next
        tcd = tc4[idxs,]
      
        x = data.frame(fam = fam, cl = cl, org = org, idr = idr, ido = ido,
          ide = tcd$id, chr = tcd$chr, beg = tcd$beg, end = tcd$end)
        dat = rbind(dat, x)
      }
    }
  }
}
write.table(dat, file = fo, sep = "\t", row.names = F, col.names = T, quote = F)
