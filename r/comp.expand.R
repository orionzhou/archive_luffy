require(ape)
require(geiger)
require(igraph)

qnames = c(
  "HM058", "HM125", "HM056", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340"
)
tname = "HM101"
chrs = sprintf("chr%s", 1:8)

dir = file.path(Sys.getenv("misc3"), "comp.ortho")

##### 
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
