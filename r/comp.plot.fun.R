require(plyr)
require(rtracklayer)
require(Cairo)
require(grid)
require(Rsamtools)

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

## data processing functions
prep_coord_mapping <- function(dcoo, seqinfo, gap_prop = 0.4, gap_len = 5000) {
  dfi = dcoo
  if(ncol(dcoo) == 3) { dfi = cbind(dcoo, srd = "+", stringsAsFactors = F) }
  colnames(dfi)[1:4] = c('chr', 'beg', 'end', 'srd')
  order.idx = order(dfi[,1], dfi[,2], dfi[,3])
  df = dfi[order.idx,]
  panels = rep(0, nrow(df))
  cnt = 1
  for (i in 1:nrow(df)) {
    if( panels[i] > 0 ) { next }
    panels[i] = cnt
    idP = df[i,1]; begP = df[i,2]; endP = df[i,3]
    for (j in (i+1):nrow(df)) {
      id = df[j,1]; beg = df[j,2]; end = df[j,3]
      if( j == i | j > nrow(df) | id != idP ) { break }
      if( end < begP ) {
        gap = begP - end - 1
      } else if( endP < beg ) {
        gap = beg - endP + 1
      } else {
        gap = 0
      }
      panelBeg = min(beg, begP)
      panelEnd = max(end, endP)
      panelLen = panelEnd - panelBeg + 1
      if( gap/panelLen <= gap_prop | gap <= gap_len ) {
        panels[j] = cnt
        idP = id; begP = panelBeg; endP = panelEnd
      } else {
        break
      }
    }
    cnt = cnt + 1
  }
  
  dcon = cbind(dfi, panel_old = panels[order(order.idx)])
  panel_unique = unique(dcon$panel_old)
  panel_mapping = data.frame(panel_old = panel_unique, 
    pan = 1:length(panel_unique))
  dcon = cbind(dcon, idx = 1:nrow(dcon))
  dcon = merge(dcon, panel_mapping, by = 'panel_old')
  dcon = dcon[order(dcon$idx),]
  dcon = dcon[, !colnames(dcon) %in% c('idx','panel_old')]
  
  tmp1 = ddply(dcon, .(pan, srd), summarise, len = sum(end - beg));
  tmp2 = ddply(tmp1, .(pan), summarise, srd = srd[which(len == max(len))][1])
  dpan = ddply(dcon, .(pan), summarise, chr = unique(chr), beg = min(beg), 
    end = max(end))
  dpan = merge(dpan, tmp2, by = 'pan')

  dlen = data.frame(id = as.character(seqnames(seqinfo)), 
    size = seqlengths(seqinfo), stringsAsFactors = F)
  dpan = merge(dpan, dlen, by.x = 'chr', by.y = 'id')
  dpan = dpan[order(dpan$pan),]
  begr = as.integer( dpan$beg - 0.02 * (dpan$end - dpan$beg) )
  endr = as.integer( dpan$end + 0.02 * (dpan$end - dpan$beg) )
  dpan = cbind(dpan, begr = begr, endr = endr)
  dpan$begr = apply(dpan, 1, function(x) max(1, as.numeric(x['begr'])))
  dpan$endr = apply(dpan, 1, function(x) 
    min(as.numeric(x['size']), as.numeric(x['endr'])))
  dpan = cbind(dpan, len = dpan$endr - dpan$begr + 1)

  len_total = sum(dpan$len)
  itv = as.integer(0.01 * len_total)
  dpan = cbind(dpan, beg.a = itv + 1)
  if(nrow(dpan) > 1) {
    for (i in 2:nrow(dpan)) {
      dpan$beg.a[i] = dpan$beg.a[i - 1] + dpan$len[i - 1] + itv + 1
    }
  }
  
  dmap = dpan[,c('pan','chr','begr','endr','srd','len','beg.a')]
  colnames(dmap)[3:4] = c('beg','end')
  cbind(dmap, end.a = dmap$beg.a + dmap$len - 1)
}
prep_ticks <- function(dmap, tick_itv) {
  dtik = data.frame(pan = c(), pos = c())  
  for (i in 1:nrow(dmap)) {
    tick_beg = pretty(c(dmap$beg[i], dmap$end[i]))[1]
    ticks = seq(tick_beg, dmap$end[i], by = tick_itv)
    ticks = ticks[ticks >= dmap$beg[i] & ticks <= dmap$end[i]]
    if(length(ticks) == 0) { ticks = tick_beg }
    dtik = rbind(dtik, data.frame(pan = rep(dmap$pan[i], length(ticks)), 
      pos = ticks))
  }
  
  dtik = merge(dtik, dmap, by = 'pan')
  dtik = cbind(dtik, pos.a = 0)
  for (i in 1:nrow(dtik)) {
    if(dtik$srd[i] == "-") {
      dtik$pos.a[i] = dtik$end.a[i] - (dtik$pos[i] - dtik$beg[i])
    } else {
      dtik$pos.a[i] = dtik$beg.a[i] + (dtik$pos[i] - dtik$beg[i])
    }
  }
  dtik = dtik[dtik$pos.a > 0, c('pan','pos','srd','pos.a')]
  dtik
}

coord_mapping <- function(dcoo, dmap) {
  if(is.null(dcoo)) { return(NULL) }
  if(ncol(dcoo) == 3) { dcoo = cbind(dcoo, srd = "+", stringsAsFactors = F) }
  colnames(dcoo)[1:4] = c('chr','beg','end','srd')
  idxs.raw = c()
  pan.idxs = c()
  for (i in 1:nrow(dmap)) {
    chr = dmap$chr[i]; beg = dmap$beg[i]; end = dmap$end[i]; pan = dmap$pan[i]
    idxss = which( dcoo$chr == chr & ( (beg <= dcoo$beg & dcoo$beg <= end) | 
      (beg <= dcoo$end & dcoo$end <= end) | 
      (beg > dcoo$beg & end < dcoo$end) ) )
    if(length(idxss) > 0) {
      idxs.raw = c(idxs.raw, idxss)
      pan.idxs = c(pan.idxs, rep(i, length(idxss)))
    }
  }
  
  if(length(idxs.raw) == 0) {idxs=NULL} else {idxs=sort(unique(idxs.raw))}
  pans = c(); begs.a = c(); ends.a = c(); srds.a = c()
  for (i in idxs) {
    chr = dcoo[i,1]; beg = dcoo[i,2]; end = dcoo[i,3]; srd = dcoo[i,4]
    
    pan.idxss = pan.idxs[ which(idxs.raw == i) ]
    pan.idx = pan.idxss[1]
    if( length(pan.idxss) > 1 ) {
      lens_ovlp = c()
      for (i in 1:length(pan.idxss) ) {
        len_ovlp = min(dmap$end[pan.idxss[i]], end) - 
          max(dmap$beg[pan.idxss[i]], beg) + 1
        lens_ovlp = c(lens_ovlp, len_ovlp)
      }
      pan.idx = pan.idxss[ which(lens_ovlp == max(lens_ovlp))[1] ]
    }
    
    pan = dmap$pan[pan.idx]
    pan.beg = dmap$beg[pan.idx]
    pan.end = dmap$end[pan.idx]
    pan.srd = dmap$srd[pan.idx]
    pan.beg.a = dmap$beg.a[pan.idx]
    pan.end.a = dmap$end.a[pan.idx]
    
    beg = max(beg, pan.beg)
    end = min(pan.end, end)
    beg.a = ifelse( pan.srd == "-", pan.end.a - (end - pan.beg), 
      pan.beg.a + (beg - pan.beg) );
    end.a = ifelse( pan.srd == "-", pan.end.a - (beg - pan.beg), 
      pan.beg.a + (end - pan.beg) );
    srd.a = ifelse(pan.srd == srd, "+", "-")
    
    pans = c(pans, pan)
    begs.a = c(begs.a, beg.a)
    ends.a = c(ends.a, end.a)
    srds.a = c(srds.a, srd.a)
  }
  cbind(dcoo[idxs,], pan = pans, beg.a = begs.a, end.a = ends.a, 
    srd.a = srds.a, stringsAsFactors = F)
}

coord_mapping_wig <- function(bw, dmap) {
  h = hash()
  for (i in 1:nrow(dmap)) {  
    chr = dmap$chr[i]; beg = dmap$beg[i]; end = dmap$end[i]; pan = dmap$pan[i]
    idxs = which( seqnames(bw) == chr & ( (beg <= start(bw) & start(bw) <= end) 
      | (beg <= end(bw) & end(bw) <= end) ) )
    h[idxs] = i
  }
  
  idxs = as.numeric(keys(h)); panel.idxs = values(h)
  bws = bw[idxs, ]
  dfi = data.frame(chr = seqnames(bws), beg = start(bws), end = end(bws), 
    srd = "+", score = score(bws), stringsAsFactors = F)
  
  panels = c(); begs.a = c(); ends.a = c(); strands.a = c()
  for (i in 1:nrow(dfi)) {
    id=dfi[i,1]; beg=dfi[i,2]; end=dfi[i,3]; strand=dfi[i,4]
    panel.idx = panel.idxs[i]
    
    panel = dmap$panel[panel.idx]
    panel.beg = dmap$beg[panel.idx]
    panel.end = dmap$end[panel.idx]
    panel.strand = dmap$strand[panel.idx]
    panel.beg.a = dmap$beg.a[panel.idx]
    panel.end.a = dmap$end.a[panel.idx]
    
    beg = max(beg, panel.beg)
    end = min(panel.end, end)
    beg.a = ifelse( panel.strand == "-", panel.end.a - (end - panel.beg), panel.beg.a + (beg - panel.beg) )
    end.a = ifelse( panel.strand == "-", panel.end.a - (beg - panel.beg), panel.beg.a + (end - panel.beg) )
    strand.a = panel.strand
    
    panels = c(panels, panel)
    begs.a = c(begs.a, beg.a)
    ends.a = c(ends.a, end.a)
    strands.a = c(strands.a, strand.a)
  }
  cbind(dfi, panel=panels, beg.a=begs.a, end.a=ends.a, stringsAsFactors=F)
}

granges2df <- function(gr) {
  ds = data.frame(chr = as.character(seqnames(gr)), beg = start(gr), 
    end = end(gr), srd = as.character(strand(gr)), stringsAsFactors = F)
  ds$srd[ds$srd == "*"] = "+"
  ds
}
prep_plot_data <- function(gro, cfgs, tname, qnames) {
  tcfg = cfgs[[tname]]
  tmap = prep_coord_mapping(granges2df(gro), tcfg$seqinfo)
  gr = GRanges(seqnames = tmap$chr, ranges = IRanges(tmap$beg, end = tmap$end),
    seqinfo = tcfg$seqinfo)
  
  dats = list()
  max_len = tmap$end.a + tmap$beg.a[1] - 1
  max_pan_len = max(tmap$len)
  for (qname in qnames) {
    cfg = cfgs[[qname]]
    aln = read_gax(cfg$tgal, cfg$tgax, cfg$tsnp, gr)

    if(is.null(aln)) {
      dats[[qname]] = NULL
      next
    }
    tg = aln$tg
    tc = aln$tc
    qmap = prep_coord_mapping(tc[,6:9], cfg$seqinfo) #tg[,5:8]
    grq = GRanges(seqnames = qmap$chr, 
      ranges = IRanges(qmap$beg, end = qmap$end), seqinfo = cfg$seqinfo)
    max_len = max(max_len, qmap$end.a + qmap$beg.a[1] - 1)
    max_pan_len = max(max_pan_len, qmap$Len)
    
    dats[[qname]] = list(tg = tg, qmap = qmap, gr = grq)
  }
  
  tick_itv = diff( pretty(c(1, max_pan_len))[1:2] )
  
  tik = prep_ticks(tmap, tick_itv)
  gap = coord_mapping(read_gap(tcfg$gap, gr), tmap)
  gene = coord_mapping(read_gene(tcfg$gene, gr), tmap)
#  mapp = coord_mapping_wig(tcfg$mapp, tmap)
  dats[[tname]] = list(map = tmap, gr = gr, 
    tik = tik, gap = gap, gene = gene)

  for (qname in qnames) {
    cfg = cfgs[[qname]]
    dat = dats[[qname]]
    if(is.null(dat)) next
    tg = dat$tg; qmap = dat$qmap; grq = dat$gr
    
    tik = prep_ticks(qmap, tick_itv)
    gap = coord_mapping(read_gap(cfg$gap, grq), qmap)
    gene = coord_mapping(read_gene(cfg$gene, grq), qmap)

    tdcoo = coord_mapping(tg[,1:4], tmap)
    qdcoo = coord_mapping(tg[,5:8], qmap)
    stopifnot(rownames(tdcoo) == rownames(tg), rownames(qdcoo) == rownames(tg))
    comp = cbind(tg, 
      tbeg.a = tdcoo$beg.a, tend.a = tdcoo$end.a, tsrd.a = tdcoo$srd.a, 
      qbeg.a = qdcoo$beg.a, qend.a = qdcoo$end.a, qsrd.a = qdcoo$srd.a, 
      stringsAsFactors = F)
    
    dats[[qname]] = list(map = qmap, gr = grq, 
      tik = tik, gap = gap, gene = gene, comp = comp)
  }
  dats$max_len = max_len
  dats
}
comp.plot <- function(fn, dats, tname, qnames, width = 1000, subtitle = "") {
  max_len = dats$max_len
  
  tdat = dats[[tname]]
  tmap = tdat$map; ttik = tdat$tik; tgap = tdat$gap; tgene = tdat$gene
  
  main = sprintf("%s compare to %d accessions", toupper(tname), length(qnames))
  fillg = c('te' = 'slategray3', 'gene' = 'tan', 'nbs' = 'forestgreen', 
    'crp' = 'dodgerblue')

  trackheight = c('axis' = 30, 'gap' = 10, 'gene' = 15, 
    'taxis' = 30, 'tgap' = 10, 'tgene' = 15, 'link' = 45)
  tracks = list()
  tracks[[tname]] = c('axis')
  for (qname in qnames) {
    tracks[[qname]] = c('tgene', 'tgap', 'taxis', 'link', 'gap', 'axis', 'gene')
  }
  
  hheight = 100
  theight = sum(trackheight[tracks[[tname]]])
  qheight = sum(trackheight[tracks[[qnames[1]]]])
  height = hheight + theight + length(qnames) * qheight
  
  lwidth = 80
  
  CairoPDF(file = fn, width = width/72, height = height/72, bg = 'transparent')
  grid.newpage()
  
  ht = 0
  vhl <- viewport(x = 0, y = unit(1, 'npc'), 
    width = unit(lwidth, 'points'), height = unit(hheight, 'points'), 
    just = c('left', 'top'), name = 'headl')
  pushViewport(vhl)
#  grid.rect(gp = gpar(fill = 'pink', alpha = 0.1, lwd = 0))
  upViewport()
  vhr <- viewport(x = unit(lwidth, 'points'), 
    y = unit(1, 'npc'), 
    width = unit(1, 'npc') - unit(lwidth, 'points'), 
    height = unit(hheight, 'points'), xscale = c(1, max_len), 
    just = c('left', 'top'), name = 'headr')
  pushViewport(vhr)
#  grid.rect(gp = gpar(fill = 'grey', alpha = 0.1, lwd = 0))
  upViewport()
  
#   grid.rect(x = 0, y = 0, width = 1, 
#     height = unit(1, 'npc') - unit(hheight, 'points'), 
#     just = c('left', 'bottom'),
#     gp = gpar(fill = NA, alpha = 1, col = 'black', lwd = 0.5))
  
  ht = hheight
  vp <- viewport(x = unit(lwidth, 'points'), 
    y = unit(1, 'npc') - unit(ht, 'points'), 
    width = unit(1, 'npc') - unit(lwidth, 'points'), 
    height = unit(1, 'npc') - unit(hheight, 'points'), 
    just = c('left', 'top'), name = 'grid')
  pushViewport(vp)
  plot_grid(width - lwidth)
  upViewport()
  
  vs = list()
  vls = list()
  vrs = list()
  
  os = hheight
  for (name in c(tname, qnames)) {
    ht = sum(trackheight[tracks[[name]]])
    vl <- viewport(x = 0, y = unit(1, 'npc') - unit(os, 'points'), 
      width = unit(lwidth, 'points'), height = unit(ht, 'points'), 
      just = c('left', 'top'), name = sprintf("%s.l", name))
    pushViewport(vl)
    upViewport()
    vr <- viewport(x = unit(lwidth, 'points'), 
      y = unit(1, 'npc') - unit(os, 'points'), 
      width = unit(1, 'npc') - unit(lwidth, 'points'), 
      height = unit(ht, 'points'), xscale = c(1, max_len), 
      just = c('left', 'top'), name = sprintf("%s.r", name))
    pushViewport(vr)
    upViewport()
    vls[[name]] = vl
    vrs[[name]] = vr
    
    grid.rect(x = 0, y = unit(1, 'npc') - unit(os, 'points'), width = 1, 
      height = unit(ht, 'points'), just = c('left', 'top'),
      gp = gpar(fill = NA, alpha = 1, col = 'black', lwd = 0.5))
    os = os + ht
  }
  
  plot_title(main, subtitle, fillg, max_len, vp = vhr)
  for (name in c(tname, qnames)) {
    vl = vls[[name]]; vr = vrs[[name]]
    cht = sum(trackheight[tracks[[name]]])
    for (track in tracks[[name]]) {
      ht = trackheight[track]
      y = unit(cht - ht / 2, 'points')
      if(track == 'axis') {
        plot_legend(sprintf("%s (kb)", name), y = y, vp = vl)
        text.above = ifelse(name == tname, F, T)
        plot_axis(dats[[name]], y = y, text.above = text.above, vp = vr)
      } else if(track == 'taxis') {
        plot_legend(sprintf("%s (kb)", tname), y = y, vp = vl)
        plot_axis(dats[[tname]], y = y, text.above = T, vp = vr)
      } else if(track == "gap") {
        plot_legend(sprintf("%s gap", name), y = y, vp = vl)
        plot_gap(dats[[name]], y = y, vp = vr)
      } else if(track == 'tgap') {
        plot_legend(sprintf("%s gap", tname), y = y, vp = vl)
        plot_gap(dats[[tname]], y = y, vp = vr)
      } else if(track == "gene") {
        plot_legend(sprintf("%s gene", name), y = y, vp = vl)
        plot_gene(dats[[name]], fillg, y = y, vp = vr)
      } else if(track == 'tgene') {
        plot_legend(sprintf("%s gene", tname), y = y, vp = vl)
        plot_gene(dats[[tname]], fillg, y = y, vp = vr)
      } else if(track == "link") {
        plot_legend(sprintf("%s/%s", name, tname), y = y, vp = vl)
        plot_link(dats[[name]], y = y, height = ht, vp = vr)
      }
      cht = cht - ht
    }
  }
  dev.off()
}

## plotting functions
plot_grid <- function(wd, itv = 30, vp = NULL) {
  n = wd %/% itv
  cols = rep('lightsteelblue3', n + 1)
  cols[1] = 'red'
  xs = unit(seq(0, length.out = n + 1, by = itv), 'points')
  grid.segments( x0 = xs, x1 = xs, y0 = 0, y1 = 1, 
    gp = gpar(col = cols, alpha = 1, lwd = 0.5), vp = vp)
}
plot_legend <- function(str, x = unit(4, 'points'), y = unit(0.5, 'npc'), 
  text.rot = 0, vp = NULL) {
  grid.text( label = str, 
    x = x, y = y, just = c("left", "center"), 
    rot = text.rot, gp = gpar(cex = 0.8, fontfamily = "Mono"), 
    vp = vp)
}
plot_axis <- function(dat, y = unit(0.5, 'npc'), col.p = 'red', col.n = 'blue',
  text.above = F, text.rot = 0, vp = NULL) {
  if(is.null(dat$map)) return(NULL)
  dax = dat$map[,c('chr', 'beg.a', 'end.a', 'srd')]
  colnames(dax) = c('id','beg','end','srd')
  
  line.cols = rep(col.p, nrow(dax))
  line.cols[which(dax$srd == "-")] = col.n
  grid.segments(
    x0 = unit(dax$beg, 'native'), x1 = unit(dax$end, 'native'),
    y0 = y, y1 = y, 
    gp = gpar(col = line.cols), vp = vp)
  
  text.x = c()
  text.just = c()
  for (i in 1:nrow(dax)) {
    beg = dax[i,2]; end = dax[i,3]; srd = dax[i,4]
    if(srd == '-') {
      x0 = unit(beg, 'native') + unit(5,'points')
      x1 = unit(beg, 'native')      
      y0 = y + unit(ifelse(text.above, -3, 3), 'points')
      y1 = y
      text.x = c(text.x, end)
      text.just = c(text.just, 'right', 'center')
    } else {
      x0 = unit(end, 'native') - unit(5,'points')
      x1 = unit(end, 'native')
      y0 = y + unit(ifelse(text.above, -3, 3), 'points')
      y1 = y
      text.x = c(text.x, beg)
      text.just = c(text.just, 'left', 'center')
    }
    if(i == 1) {
      arrows.x0 = x0; arrows.x1 = x1
      arrows.y0 = y0; arrows.y1 = y1
    } else {
      arrows.x0 = unit.c(arrows.x0, x0); arrows.x1 = unit.c(arrows.x1, x1)
      arrows.y0 = unit.c(arrows.y0, y0); arrows.y1 = unit.c(arrows.y1, y1)
    }
  }
  grid.segments( 
    x0 = arrows.x0, x1 = arrows.x1, y0 = arrows.y0, y1 = arrows.y1,
    gp = gpar(col = line.cols), vp = vp)
  
  text.y = y + unit(ifelse(text.above, 5, -5), 'points') 
  text.offset = y + unit(ifelse(text.above, -10, 10), "points")
  grid.text( label = dax$id, 
    x = unit((dax$beg + dax$end) / 2, 'native'), 
    y = text.y, just = c("center", "center"), 
    rot = text.rot, gp = gpar(cex = 0.7, fontfamily = "mono"), vp = vp)
  
  dtik = dat$tik
#   dlin = ddply(dtik, .(pan), summarise, beg.a = min(pos.a), end.a = max(pos.a))
#   grid.segments(
#     x0 = unit(dlin$beg.a, 'native'),
#     x1 = unit(dlin$end.a, 'native'),
#     y0 = y, y1 = y, vp = vp)
  
  if( !text.above ) {
    tick.y = y + unit(3, 'points')
    text.y = y + unit(5, 'points')
    text.just = c('center', 'bottom')
  } else {
    tick.y = y - unit(3, 'points')
    text.y = y - unit(5, 'points')
    text.just = c('center', 'top')
  }
  grid.segments(
    x0 = unit(dtik$pos.a, 'native'),
    x1 = unit(dtik$pos.a, 'native'),
    y0 = y, y1 = tick.y, 
    vp = vp)
  grid.text( 
    label = dtik$pos / 1000,
    x = unit(dtik$pos.a, 'native'), 
    y = text.y,  just = text.just, 
    gp = gpar(cex = 0.6, fontfamily = "mono"), vp = vp)
}
plot_gene <- function(dat, fillg, y = unit(0.5, 'npc'), text.show = F, 
  text.offset = unit(10, 'points'), text.above = F, 
  text.rot = 0, vp = NULL) {
  if(is.null(dat$gene)) { return(NULL) }
  dgen = dat$gene[,c('beg.a', 'end.a', 'srd', 'type', 'cat')]
  colnames(dgen) = c('beg', 'end', 'srd', 'type', 'cat')
  
  grid.segments(x0 = unit(min(dgen$beg), 'native'), 
    x1 = unit(max(dgen$end), 'native'), y0 = y, y1 = y, 
    gp = gpar(col = 'grey', lty = 3), vp = vp)
  
  dfp = dgen[dgen$srd == "+",]
  dfn = dgen[dgen$srd == "-",]
  height = 5
  yp = y + unit(height / 2, 'points')
  yn = y - unit(height / 2, 'points')
  
  if( !empty(dfp) ) {
    for (cat in names(fillg)) {
      dfm = dfp[dfp$type == 'mrna' & dfp$cat == cat,]
      dfc = dfp[dfp$type == 'cds' & dfp$cat == cat,]
      if(!empty(dfm)) {
        grid.segments(
        x0 = unit(dfm$beg, 'native'), x1 = unit(dfm$end, 'native'),
        y0 = yp, y1 = yp, 
        gp = gpar(col = fillg[dfm$cat]), vp = vp)
      }
      if(!empty(dfc)) {
        grid.rect( 
        x = unit(dfc$beg, 'native'), y = yp,
        width = unit(dfc$end - dfc$beg, 'native'), 
        height = unit(height, 'points'), just = c('left', 'center'),
        gp = gpar(lwd = 0, fill = fillg[dfc$cat], alpha = 0.9), vp = vp)
      }
    }
  }
  
  if( !empty(dfn) ) {
    for (cat in names(fillg)) {
      dfm = dfn[dfn$type == 'mrna' & dfn$cat == cat,]
      dfc = dfn[dfn$type == 'cds' & dfn$cat == cat,]
      if(!empty(dfm)) {
        grid.segments(
        x0 = unit(dfm$beg, 'native'), x1 = unit(dfm$end, 'native'),
        y0 = yn, y1 = yn, 
        gp = gpar(col = fillg[dfm$cat]), vp = vp)
      }
      if(!empty(dfc)) {
        grid.rect( 
        x = unit(dfc$beg, 'native'), y = yn,
        width = unit(dfc$end - dfc$beg, 'native'), 
        height = unit(height, 'points'), just = c('left', 'center'),
        gp = gpar(lwd = 0, fill = fillg[dfc$cat], alpha = 0.9), vp = vp)
      }
    }
  }

  if( text.show ) {
    text.y = ifelse(text.above, y + text.offset, y - text.offset)
    grid.text( df$id, 
      x = unit(df$beg.a, 'native'), y = text.y, just = c("left", "center"), 
      rot = text.rot, gp = gpar(cex = 0.8, fontfamily = "sans"), 
      vp = vp)
  }
}
plot_gap <- function(dat, y = unit(0.5, 'npc'), fill = 'grey', vp = NULL) {
  if(is.null(dat$gap)) { return(NULL) }
  dgap = dat$gap[,c('chr', 'beg.a', 'end.a')]
  colnames(dgap) = c('id','beg','end')
  height = unit(5, 'points')
  grid.rect( 
    x = unit(dgap$beg, 'native'), y = y,
    width = unit(dgap$end - dgap$beg, 'native'), height = height,
    just = c('left', 'center'),
    gp = gpar(lwd = 0, fill = fill, alpha = 0.9), vp = vp)
}
plot_link <- function(dat, y = unit(0.5, 'npc'), height = 30, 
  fill.p = 'skyblue1', fill.n = 'tomato', alpha = 0.5, vp = NULL) {
  comp = dat$comp
  if(is.null(comp)) { return(NULL) }
  comp.xs = c()
  comp.fills = c()
  comp.ids = c()
  for (i in 1:nrow(comp) ) {
    if( comp$qsrd.a[i] == comp$tsrd.a[i] ) {
      comp.x = c(comp$tbeg.a[i], comp$tend.a[i], comp$qend.a[i], comp$qbeg.a[i])
      comp.fill = fill.p
      comp.id = rep(i, 4)
    } else {
      comp.x = c(comp$tbeg.a[i], comp$tend.a[i], comp$qbeg.a[i], comp$qend.a[i])
      comp.fill = fill.n
      comp.id = rep(i, 4)
    }
    comp.xs = c(comp.xs, comp.x)
    comp.fills = c(comp.fills, comp.fill)
    comp.ids = c(comp.ids, comp.id)
  }
  
  tmp = rep(c(1, 1, -1, -1), nrow(comp))
  grid.polygon(
    x = unit(comp.xs, 'native'),
    y = y + unit(tmp * height / 2, 'points'),
    id = comp.ids,
    gp = gpar(fill = comp.fills, alpha = alpha, lwd = 0),
    vp = vp)
}
plot_hist <- function(df, fill='grey', vp=NULL) {
  colnames(df)[1:3] = c("beg", "end", "score")
  grid.rect( 
    x = unit(df$beg,'native'), y=unit(0, 'native'), 
    width = unit(df$end-df$beg,'native'), height=unit(df$score,'native'), 
    just=c('left','bottom'),
    gp = gpar(lwd=0, fill=fill, alpha=1), vp=vp)
}
plot_title <- function(main, subtitle, fill, max_len, vp = NULL) {
  grid.text(main, 
    x = unit(0.5, 'npc'), y = unit(0.7, 'npc'), 
    gp = gpar(fontface = 'bold', fontfamily = 'serif'),
    just = c("center", "top"), 
    vp = vp)
  grid.text(subtitle, 
    x = unit(0.5, 'npc'), y = unit(0.7, 'npc') - unit(1, 'lines'), 
    gp = gpar(fontface = 'bold', fontfamily = 'serif'),
    just = c("center", "top"), 
    vp = vp)
  
  n = length(fill)
  fill.labels = names(fill)
  grid.rect( 
    x = rep(unit(0.02, 'npc'), n), 
    y = unit(0.1, 'npc') + unit(seq(0, by = 10, length.out = n), 'points'),
    width = unit(30, 'points'), height = unit(5, 'points'),
    just = c('left', 'bottom'),
    gp = gpar(lwd = 0, fill = fill, alpha = 0.9), vp = vp)
  grid.text( fill.labels, 
    x = unit(0.02, 'npc') + unit(40, 'points'), 
    y = unit(0.1, 'npc') + unit(seq(0, by = 10, length.out = n), 'points'), 
    just = c("left", "bottom"), 
    gp = gpar(cex = 0.9, fontfamily = "serif"), vp = vp)

  len = diff( pretty(1:max_len, 20)[1:2] )
  name = sprintf("%.0fkb", len/1000)
  grid.segments( 
    x0 = unit(0.98, 'npc') - unit(len, 'native'), x1 = unit(0.98, 'npc'),
    y0 = 0.2, y1 = 0.2, 
    vp = vp)
  grid.segments( 
    x0 = unit.c(unit(0.98, 'npc') - unit(len, 'native'), unit(0.98, 'npc')),
    x1 = unit.c(unit(0.98, 'npc') - unit(len, 'native'), unit(0.98, 'npc')),
    y0 = rep(0.2, 2), 
    y1 = rep(unit(0.2, 'npc') + unit(3, 'points'), 2),
    vp = vp)
  grid.text( name,
    x = unit(0.98, 'npc') - unit(len / 2, 'native'), 
    y = unit(0.2, 'npc') + unit(5, 'points'), 
    just = c("center", "bottom"), 
    gp = gpar(cex = 0.9, fontfamily = "serif"), vp = vp)
}


