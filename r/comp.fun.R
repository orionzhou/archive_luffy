require(plyr)
require(rtracklayer)
require(Rsamtools)

get_genome_cfg <- function(org) {
  gdir = file.path(Sys.getenv("genome"), toupper(org))
  
  size = file.path(gdir, "15.sizes")
  
  gap = file.path(gdir, "16.gap.bed")
  gapz = file.path(gdir, "16.gap.bed.gz")
  
  gene = file.path(gdir, "51.tbl")
  genez = file.path(gdir, "51.tbl.gz")

  dna = file.path(gdir, '11_genome.fas')
  protein = file.path(gdir, '51.fas')

  list(size = size, gap = gap, gapz = gapz, gene = gene, genez = genez,
    dna = dna, protein = protein)
}
read_seqinfo <- function(fsize) {
  tsize = read.table(fsize, header = F, sep = "\t",  as.is = T, 
    col.names = c("id", "size"))
  Seqinfo(tsize$id, seqlengths = tsize$size)
}
get_comp_cfg <- function(tname, qnames) {
  cfgs = list()
  
  gdir = file.path(Sys.getenv("genome"), toupper(tname))
  seqinfo = read_seqinfo(file.path(gdir, "15.sizes"))
  gap = file.path(gdir, "16.gap.bed.gz")
  gene = file.path(gdir, "51.tbl.gz")
  cfgs[[tname]] = list(seqinfo = seqinfo, gap = gap, gene = gene)
  
  for (qname in qnames) {
    gdir = file.path(Sys.getenv("genome"), toupper(qname))
    seqinfo = read_seqinfo(file.path(gdir, "15.sizes"))
    gap = file.path(gdir, "16.gap.bed.gz")
    gene = file.path(gdir, "51.tbl.gz")

    cdir = sprintf("%s/%s_%s/23_blat", Sys.getenv("misc3"), 
      toupper(qname), toupper(tname))
    tgal = sprintf("%s/31.9/gal.gz", cdir)
    qgal = sprintf("%s/41.9/gal.gz", cdir)
    tgax = sprintf("%s/31.9/gax.gz", cdir)
    qgax = sprintf("%s/41.9/gax.gz", cdir)
    tsnp = sprintf("%s/31.9/snp.gz", cdir)
    qsnp = sprintf("%s/41.9/snp.gz", cdir)
    
    vdir = sprintf("%s/hapmap_40/30_vnt/%s", Sys.getenv("misc3"), toupper(qname))
    fsnp = sprintf("%s/snp.gz", vdir)
    fidm = sprintf("%s/idm.gz", vdir)
    fhet = sprintf("%s/het.gz", vdir)
    fcov = sprintf("%s/../../11_pipe_mapping/35_cov/%s.bw", vdir, toupper(qname))
    fcovab = sprintf("%s/../../40_sv/01_ab/%s.bw", vdir, toupper(qname))
    
    cfgs[[qname]] = list(seqinfo = seqinfo, gap = gap, gene = gene,
      cdir = cdir, 
      tgal = tgal, qgal = qgal, tgax = tgax, tsnp = tsnp, qsnp = qsnp,
      vdir = vdir)
  }
  cfgs
}

parse_tabix <- function(txt) 
  read.csv(textConnection(txt), sep = "\t", header = F, stringsAsFactors = F)
trim_gax <- function(ds, beg, end) {
  for (i in 1:nrow(ds)) {
    if(ds$tbeg[i] < beg) {
      if(as.character(ds$tsrd[i]) == as.character(ds$qsrd[i])) {
        ds$qbeg[i] = ds$qbeg[i] + (beg - ds$tbeg[i])
      } else {
        ds$qend[i] = ds$qend[i] - (beg - ds$tbeg[i])
      }
      ds$tbeg[i] = beg
    }
    if(end < ds$tend[i]) {
      if(as.character(ds$tsrd[i]) == as.character(ds$qsrd[i])) {
        ds$qend[i] = ds$qend[i] - (ds$tend[i] - end)
      } else {
        ds$qbeg[i] = ds$qbeg[i] + (ds$tend[i] - end)
      }
      ds$tend[i] = end
    }
  }
  ds
}
read_gax <- function(fgal, fgax, fsnp, gr) {
  gr = reduce(gr)
  
  tg = data.frame()
  tc = data.frame()
  
  gax = open(TabixFile(fgax))
  grs = gr[as.character(seqnames(gr)) %in% seqnamesTabix(gax)]
  if(length(grs) == 0) return(NULL)
  x = scanTabix(gax, param = grs)
  close(gax)
  
  txts = rapply(x, c)
  if(length(txts) == 0) return(NULL)
  
  lc = list()
  for (i in 1:length(grs)) {
    if(length(x[[i]]) == 0) next
    df1 = parse_tabix(x[[i]])
    colnames(df1) = c('tid', 'tbeg', 'tend', 'tsrd', 
      'qid', 'qbeg', 'qend', 'qsrd', 'cid', 'lev')
    tgs = trim_gax(df1, start(grs)[i], end(grs)[i])
    tgs = cbind(tgs, len = tgs$qend - tgs$qbeg + 1)
    
    tcs = ddply(tgs, .(cid), summarise, 
      tid = unique(tid), tbeg = min(tbeg), tend = max(tend), 
      tsrd = unique(tsrd),
      qid = unique(qid), qbeg = min(qbeg), qend = max(qend), 
      qsrd = unique(qsrd),
      ali = sum(len), gapo = length(len) - 1)
    for (cid in tcs$cid) {
      if(is.null(lc[[as.character(cid)]])) {
        lc[[as.character(cid)]] = 0
        next
      }
      lc[[as.character(cid)]] = lc[[as.character(cid)]] + 1
      ncid = cid + lc[[as.character(cid)]] / 100
      tcs$cid[tcs$cid == cid] = ncid
      tgs$cid[tgs$cid == cid] = ncid
    }
    tg = rbind(tg, tgs)
    tc = rbind(tc, tcs)
  }
  
  tc = cbind(tc, mis = 0)
  snp = open(TabixFile(fsnp))
  for (i in 1:nrow(tc)) {
    tid = as.character(tc$tid[i]); tbeg = tc$tbeg[i]; tend = tc$tend[i]
    if(! tid %in% seqnamesTabix(snp)) next
    grs = GRanges(seqname = tid, ranges = IRanges(tbeg, end = tend))
    x = scanTabix(snp, param = grs)[[1]]
    if(length(x) == 0) next
    tso = parse_tabix(x)
    colnames(tso) = c('tid', 'tpos', 'ref', 'alt', 'qid', 'qpos', 'cid', 'lev')
    tc$mis[i] = nrow(tso)
  }
  close(snp)

#   gal = open(TabixFile(fgal))
#   grs = gr[as.character(seqnames(gr)) %in% seqnamesTabix(gal)]
#   x = scanTabix(gal, param = grs)
#   tco = parse_tabix(rapply(x, c))[,1:19]
#   colnames(tco) = c('cid', 'tid', 'tbeg', 'tend', 'tsrd', 'tsize',
#       'qid', 'qbeg', 'qend', 'qsrd', 'qsize', 
#       'lev', 'ali', 'mat', 'mis', 'qN', 'tN', 'ident', 'score')
#   close(gal)

  qgap = tc$qend - tc$qbeg + 1 - tc$ali
  tgap = tc$tend - tc$tbeg + 1 - tc$ali
  score = tc$ali * 1 + tc$mis * (-2) + tc$gapo * (-5) + 
    (qgap + tgap - tc$gapo) * (-2)
  tc = cbind(tc, score = score)
  tc = tc[order(tc$tid, tc$tbeg, tc$tend),]
  
  idxs = tc$ali >= sum(tc$ali) / 200
  cids = tc$cid[idxs]
  list(tg = tg[tg$cid %in% cids,], tc = tc[idxs,])
}
read_gap <- function(ftbx, gr) {
  gr = reduce(gr)
  tbx = open(TabixFile(ftbx))
  grs = gr[as.character(seqnames(gr)) %in% seqnamesTabix(tbx)]
  if(length(grs) == 0) return(NULL)
  x = scanTabix(tbx, param = grs)
  txts = rapply(x, c)
  close(tbx)
  if(length(txts) == 0) return(NULL)
  res = parse_tabix(txts)
  colnames(res) = c('chr', 'beg', 'end')
  res
}
read_gene <- function(ftbx, gr) {
  gr = reduce(gr)
  tbx = open(TabixFile(ftbx))
  grs = gr[as.character(seqnames(gr)) %in% seqnamesTabix(tbx)]
  if(length(grs) == 0) return(NULL)
  x = scanTabix(tbx, param = grs)
  txts = rapply(x, c)
  close(tbx)
  if(length(txts) == 0) return(NULL)
  res = parse_tabix(txts)
  colnames(res) = c('chr', 'beg', 'end', 'srd', 'id', 'type', 'cat')
  res
}


