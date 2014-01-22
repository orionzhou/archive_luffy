library(plyr)
library(rtracklayer)
library(Cairo)

read_genome_stat <- function(name) {
  dir = file.path("/home/youngn/zhoup/Data/genome", toupper(name))
  
  # sequence lengths
  f_len = file.path(dir, "15.sizes")
  len = read.table(f_len, header=F, sep="\t", as.is=T)
  colnames(len) = c("id", "size")

  # gap locations
  f_gap = file.path(dir, "16_gap.tbl")
  gap = read.table(f_gap, header=T, sep="\t", as.is=T)

  # gene annotation
  f_gen = file.path(dir, "51.tbl")
  gene = read.table(f_gen, header=T, sep="\t", quote="", as.is=T)

  list( name=name, dir=dir, len=len, gap=gap, gene=gene )
}

build_ideogram <- function(tgap, tlen) {
  gg = GRanges(seqnames=Rle(tgap$id), ranges=IRanges(tgap$beg, end=tgap$end))
  gt = GRanges(seqnames=Rle(tlen$id), ranges=IRanges(1, end=tlen$size))
  gb = setdiff(gt, gg)

  tb = data.frame(id=seqnames(gb), beg=start(gb), end=end(gb), col='gpos50')
  tg = cbind(tgap[,-4], col='gneg')
  tid = rbind(tb, tg)
  tid$beg = tid$beg-1
  tid = tid[order(tid$id, tid$beg),]
  colnames(tid) = c('chrom', 'chromStart', 'chromEnd', 'gieStain')
  cbind(tid,name=sprintf("band%d",1:nrow(tid)))
}

read_comp_stat <- function(qname, tname) {
  dir = sprintf("/home/youngn/zhoup/Data/misc3/%s_%s", toupper(qname), 
    toupper(tname))
  tw = read.table(file.path(dir, "23_blat/26.gal"), header=TRUE, sep="\t", 
    as.is=T)[,1:17]
  tl = read.table(file.path(dir, "23_blat/26.gall"), header=TRUE, sep="\t", 
    as.is=T)
  list(dir=dir, tw=tw, tl=tl)
}

get_ins <- function(t) {  
  t = t[order(t$tId, t$tBeg),]
  if(nrow(t) == 1) {
    data.frame()
  } else {
    t1 = t[-nrow(t),]
    t2 = t[-1,]
    tIns = t2$tBeg - t1$tEnd - 1
    if(t1$qSrd[1] != t1$tSrd[1]) {
      qIns = t1$qEnd-t2$qBeg-1
    } else {
      qIns = t2$qBeg-t1$qEnd-1
    }
    data.frame(id=t1$id, qId=t1$qId, tId=t1$tId, 
      tPos=ceiling(t1$tEnd+t2$tBeg)/2, 
      tIns=tIns, qIns=qIns, txt=sprintf("%d^%d", tIns, qIns))
  }
}

comp_plot <- function(chr, beg, end, fo, q, t, c, tid, f_mapp) {
  axisTrack <- GenomeAxisTrack(cex=1.0, exponent=3)
  ideoTrack <- IdeogramTrack(genome = t$name, chromosome = chr, bands = tid, 
    bevel = 0, showId = T, cex = 1.0)
  mappTrack <- DataTrack(range=f_mapp, genome=t$name, type='horizon', 
    chromosome = chr, name='mapp', background.title="tan1")
  grTrack <- GeneRegionTrack(t$gene, genome = t$name, chromosome = chr, 
    name = "genes", showId=T, cex=1.0, background.title="tan1")

  tw = c$tw
  tl = c$tl
  #tws = tw[tw$tId==chr & ( beg<=tw$tEnd & tw$tBeg<=end ), ]
  #tls = tl[tl$id %in% tws$id,]

  tls = tl[tl$tId==chr & ( beg<=tl$tEnd & tl$tBeg<=end ), ]

  if(empty(tls)) {
    cTrack <- GeneRegionTrack(genome=t$name, chromosome=chr, name=q$name, 
      showId=F, fill='dodgerblue', background.title="brown")
  } else {
    tx = data.frame(chromosome=tls$tId, start=tls$tBeg, end=tls$tEnd, 
      width=tls$tEnd-tls$tBeg+1, strand=tls$tSrd, 
      feature='protein_coding', gene=tls$id, 
      exon=sprintf("%s.%d", tls$id, 1:nrow(tls)), 
      transcript=tls$id, symbol=tls$qId)
    cTrack <- GeneRegionTrack(tx, genome=t$name, chromosome=chr, 
      name=q$name, howId=F, fill='dodgerblue', background.title="brown")
  }

  tins = ddply(tls, .(id), get_ins)

  tis = tins[tins$qIns<100 & tins$tIns<100,]
  til = tins[tins$qIns>=100 | tins$tIns>=100,]
  if(empty(tins) | empty(tis)) {
    siTrack <- GeneRegionTrack(genome=t$name, chromosome=chr, 
      name="indel-s", showId=F, showExonId=T, fill='lightgreen', 
      fontcolor='black', background.title="brown", cex=0.8, rotation=90)
  } else {
    tisz = data.frame(chromosome=tis$tId, start=tis$tPos, end=tis$tPos, 
      feature='protein_coding', gene=tis$id, exon=tis$txt, 
      transcript=tis$id, symbol=tis$id)
    siTrack <- GeneRegionTrack(tisz, genome=t$name, chromosome=chr, 
      name="indel-s", showId=F, showExonId=T, 
      fill='lightgreen', fontcolor='black', background.title="brown", 
      cex=0.8, rotation=90)
  }
  if(empty(tins) | empty(til)) {
    liTrack <- GeneRegionTrack(genome=t$name, chromosome=chr, 
      name="indel-l", showId=F, showExonId=T, 
      fill='lightseagreen', fontcolor='black', 
      background.title="brown", cex=0.8, rotation=0)
  } else {
    tilz = data.frame(chromosome=til$tId, start=til$tPos, end=til$tPos, 
      feature='protein_coding', gene=til$id, exon=til$txt, transcript=til$id, 
      symbol=til$id)
    liTrack <- GeneRegionTrack(tilz, genome=t$name, chromosome=chr, 
      name="indel-l", showId=F, showExonId=T, 
      fill='lightseagreen', fontcolor='black', background.title="brown", 
      cex=0.8, rotation=0)
  }
  CairoPNG(filename = fo, width = 1000, height = 400)
  plotTracks(
    list(ideoTrack, axisTrack, mappTrack, grTrack, cTrack, 
      siTrack, liTrack), 
    extend.left=(end-beg)/20, extend.right=(end-beg)/20, 
    from=beg, to=end, sizes=c(3,3,2,3,3,2,2))
  dev.off()
}

##### OUTDATED STUFF
## data processing functions
assign_panel <- function(dfi, gap_prop=0.4, gap_len=5000) {
  if(ncol(dfi) == 3) { dfi = cbind(dfi, strand="+") }
  colnames(dfi)[1:4] = c('id', 'beg', 'end', 'strand')
  order.idx = order(dfi[,1], dfi[,2], dfi[,3])
  df = dfi[order.idx,]
  panels = rep(0, nrow(df))
  cnt = 1
  for (i in 1:nrow(df)) {
    if( panels[i] > 0 ) { next }
    panels[i] = cnt
    idP=df[i,1]; begP=df[i,2]; endP=df[i,3]
    for (j in (i+1):nrow(df)) {
      id=df[j,1]; beg=df[j,2]; end=df[j,3]
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
        idP=id; begP=panelBeg; endP=panelEnd
      } else {
        break
      }
    }
    cnt = cnt + 1
  }
  dfo = cbind(dfi, panel_old=panels[order(order.idx)])
  panel_unique = unique(dfo$panel_old)
  panel_mapping = data.frame(panel_old=panel_unique, 
    panel=1:length(panel_unique))
  dfo = cbind(dfo, idx=1:nrow(dfo))
  dfo = merge(dfo, panel_mapping, by='panel_old')
  dfo = dfo[order(dfo$idx),]
  dfo = dfo[, !colnames(dfo) %in% c('idx','panel_old')]
  ddply(dfo, .(panel), summarise, id=unique(id), beg=min(beg), end=max(end), strand=strand[which(end-beg == max(end-beg))[1]])
}

assign_block_mapping <- function(dfi, gap_len_q=10000, gap_len_h=100000) {
  order.idx = order(dfi[,1], dfi[,2], dfi[,3])
  df = dfi[order.idx,]
  blocks = rep(0, nrow(df))
  cnt = 1
  for (i in 1:nrow(df)) {
    if( blocks[i] > 0 ) { next }
    blocks[i] = cnt
    qIdP=df$qId[i]; qBegP=df$qBeg[i]; qEndP=df$qEnd[i]; strdP=df$strand[i];
    hIdP=df$hId[i]; hBegP=df$hBeg[i]; hEndP=df$hEnd[i]
    for (j in (i+1):nrow(df)) {
      qId=df$qId[j]; qBeg=df$qBeg[j]; qEnd=df$qEnd[j]; strd=df$strand[j];
      hId=df$hId[j]; hBeg=df$hBeg[j]; hEnd=df$hEnd[j]
      if( j == i | j > nrow(df) ) { next }
      if( qId != qIdP ) {
        break
      } else if( abs(qBeg-qEndP)<gap_len_q & strd==strdP & (hId==hIdP) & ( (strd=="+" & abs(hBeg-hEndP)<gap_len_h) | (strd=="-" & abs(hEnd-hBegP)<gap_len_h) ) ) {
        blocks[j] = cnt
        qIdP=qId; qBegP=qBeg; qEndP=qEnd; strdP=strd; hIdP=hId; hBegP=hBeg; hEndP=hEnd;
      }
    }
    cnt = cnt + 1
  }
  dfo = cbind(dfi, block=blocks[order(order.idx)])
  dfb = ddply(dfo, .(block), summarise, 
    qId = unique(qId), qBeg = min(qBeg), qEnd = max(qEnd), strand = unique(strand),
    hId = unique(hId), hBeg = min(hBeg), hEnd = max(hEnd), qLen = max(qEnd)-min(qBeg)+1,
    hLen = max(hEnd)-min(hBeg)+1, qLen_aln = sum(qLen), hLen_aln = sum(hLen),
    pct = mean(pct), e = min(e), score = sum(score))
  list('df'=dfo, 'dfb'=dfb)
}

prepare_coord <- function(dfi, tl) {
  dfi = merge(dfi, tl, by='id')
  df = dfi[order(dfi$panel),]
  begr = as.integer( df$beg - 0.02 * (df$end - df$beg) )
  endr = as.integer( df$end + 0.02 * (df$end - df$beg) )
  df2 = cbind(df, begr=begr, endr=endr)
  df2$begr = apply(df2, 1, function(x) max(1,as.numeric(x['begr'])))
  df2$endr = apply(df2, 1, function(x) min(as.numeric(x['size']),as.numeric(x['endr'])))
  df3 = cbind(df2, len=df2$endr-df2$begr+1)

  len_total = sum(df3$len)
  itv = as.integer(0.01 * len_total)
  df5 = cbind(df3, beg.a = itv+1)
  if(nrow(df5) > 1) {
    for (i in 2:nrow(df5)) {
      df5$beg.a[i] = df5$beg.a[i-1] + df5$len[i-1] + itv + 1
    }
  }
  df6 = df5[,c('panel','id','begr','endr','strand','len','beg.a')]
  colnames(df6)[3:4] = c('beg','end')
  df7 = cbind(df6, end.a=df6$beg.a+df6$len-1)
  df7
}

filter_assign_panel <- function(dfi, dfb) {
  if(ncol(dfi) == 3) { dfi = cbind(dfi, strand="+") }
  colnames(dfi)[1:4] = c('id','beg','end','strand')
  idxs.raw = c()
  panel.idxs = c()
  for (i in 1:nrow(dfb)) {
    id=dfb$id[i]; beg=dfb$beg[i]; end=dfb$end[i]; blk=dfb$panel[i]
    idxss = which( dfi$id==id & ( (beg<=dfi$beg & dfi$beg<=end) | (beg<=dfi$end & dfi$end<=end) ) )
    if(length(idxss) > 0) {
      idxs.raw = c(idxs.raw, idxss)
      panel.idxs = c(panel.idxs, rep(i, length(idxss)))
    }
  }
  
  if(length(idxs.raw) == 0) {idxs=NULL} else {idxs=sort(unique(idxs.raw))}
  panels = c()
  begs.a = c()
  ends.a = c()
  strands.a = c()
  for (i in idxs) {
    id=dfi[i,1]; beg=dfi[i,2]; end=dfi[i,3]; strand=dfi[i,4]
    
    panel.idxss = panel.idxs[ which(idxs.raw == i) ]
    panel.idx = panel.idxss[1]
    if( length(panel.idxss) > 1 ) {
      lens_ovlp = c()
      for (i in 1:length(panel.idxss) ) {
        len_ovlp = min(dfb$end[panel.idxss[i]], end) - max(dfb$beg[panel.idxss[i]], beg) + 1
        lens_ovlp = c(lens_ovlp, len_ovlp)
      }
      panel.idx = panel.idxss[ which(lens_ovlp == max(lens_ovlp))[1] ]
    }
    
    panel = dfb$panel[panel.idx]
    panel.beg = dfb$beg[panel.idx]
    panel.end = dfb$end[panel.idx]
    panel.strand = dfb$strand[panel.idx]
    panel.beg.a = dfb$beg.a[panel.idx]
    panel.end.a = dfb$end.a[panel.idx]
    
    beg = max(beg, panel.beg)
    end = min(panel.end, end)
    beg.a = ifelse( panel.strand == "-", panel.end.a - (end - panel.beg), panel.beg.a + (beg - panel.beg) );
    end.a = ifelse( panel.strand == "-", panel.end.a - (beg - panel.beg), panel.beg.a + (end - panel.beg) );
    strand.a = ifelse(panel.strand == strand, "+", "-")
    
    panels = c(panels, panel)
    begs.a = c(begs.a, beg.a)
    ends.a = c(ends.a, end.a)
    strands.a = c(strands.a, strand.a)
  }
  cbind(dfi[idxs,], panel=panels, beg.a=begs.a, end.a=ends.a, strand.a=strands.a, stringsAsFactors=F)
}

filter_assign_panel_wig <- function(bw, dfb) {
  h = hash()
  for (i in 1:nrow(dfb)) {  
    id=dfb$id[i]; beg=dfb$beg[i]; end=dfb$end[i]; panel=dfb$panel[i]
    idxs = which( seqnames(bw)==id & ( (beg<=start(bw) & start(bw)<=end) | (beg<=end(bw) & end(bw)<=end) ) )
    h[idxs] = i
  }
  
  idxs = as.numeric(keys(h)); panel.idxs = values(h)
  bws = bw[idxs, ]
  dfi = data.frame(id=seqnames(bws), beg=start(bws), end=end(bws), strand="+", score=score(bws), stringsAsFactors=F)
  
  panels=c(); begs.a = c(); ends.a = c(); strands.a = c()
  for (i in 1:nrow(dfi)) {
    id=dfi[i,1]; beg=dfi[i,2]; end=dfi[i,3]; strand=dfi[i,4]
    panel.idx = panel.idxs[i]
    
    panel = dfb$panel[panel.idx]
    panel.beg = dfb$beg[panel.idx]
    panel.end = dfb$end[panel.idx]
    panel.strand = dfb$strand[panel.idx]
    panel.beg.a = dfb$beg.a[panel.idx]
    panel.end.a = dfb$end.a[panel.idx]
    
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

get_ticks <- function(dfi, tick_itv) {
  dft = data.frame(panel=c(),pos=c())  
  for (i in 1:nrow(dfi)) {
    tick_beg = pretty(c(dfi$beg[i], dfi$end[i]))[1]
    ticks = seq(tick_beg, dfi$end[i], by=tick_itv)
    if(length(ticks)==0) { ticks = tick_beg }
    dft = rbind(dft, data.frame(panel=rep(dfi$panel[i], length(ticks)), pos=ticks))
  }
  
  dft = merge(dft, dfi, by='panel')
  dft = cbind(dft, pos.a=0)
  for (i in 1:nrow(dft)) {
    if(dft$strand[i] == "-") {
      dft$pos.a[i] = dft$end.a[i] - (dft$pos[i] - dft$beg[i])
    } else {
      dft$pos.a[i] = dft$beg.a[i] + (dft$pos[i] - dft$beg[i])
    }
  }
  dft = dft[dft$pos.a > 0, c('panel','pos','strand','pos.a')]
  dfl = ddply(dft, .(panel), summarise, beg.a=min(pos.a), end.a=max(pos.a))
  list(tick=dft, line=dfl)
}


data_preprocess <- function(tws, tls, q, t) {
  qPanel = prepare_coord( assign_panel(tws[,c('qId', 'qBeg', 'qEnd', 'qSrd')]), q$len )
  tPanel = prepare_coord( assign_panel(tws[,c('tId', 'tBeg', 'tEnd', 'tSrd')]), t$len )
  
  max_len = max(qPanel$end.a+qPanel$beg.a[1]-1, tPanel$end.a+tPanel$beg.a[1]-1)
  tick_itv = diff( pretty(c(1,max(qPanel$len, tPanel$len)))[1:2] )
  xaxis1 = get_ticks(qPanel, tick_itv)
  xaxis2 = get_ticks(tPanel, tick_itv)

  dfm.a = filter_assign_panel(tls[,c('qId','qBeg','qEnd','qSrd')], qPanel)
  dfm.b = filter_assign_panel(tls[,c('tId','tBeg','tEnd','tSrd')], tPanel)
  stopifnot(rownames(dfm.a) == rownames(tls), rownames(dfm.b) == rownames(tls))
  comp = cbind(tls, qPanel=dfm.a$panel, qPanel.beg=dfm.a$beg.a, qPanel.end=dfm.a$end.a, qPanel.strand=dfm.a$strand.a, tPanel=dfm.b$panel, tPanel.beg=dfm.b$beg.a, tPanel.end=dfm.b$end.a, tPanel.strand=dfm.b$strand.a, stringsAsFactors=F)
  
  if(empty(q$gap)) {gap1=NULL} else {gap1=filter_assign_panel(q$gap, qPanel)}
  if(empty(t$gap)) {gap2=NULL} else {gap2=filter_assign_panel(t$gap, tPanel)}

  if(empty(q$gene)) {gene1=NULL} else {gene1=filter_assign_panel(q$gene, qPanel)}
  if(empty(t$gene)) {gene2=NULL} else {gene2=filter_assign_panel(t$gene, tPanel)}
  
  if(empty(q$mapp)) {mapp1=NULL} else {mapp1=filter_assign_panel_wig(q$mapp, qPanel)}
  if(empty(t$mapp)) {mapp2=NULL} else {mapp2=filter_assign_panel_wig(t$mapp, tPanel)}
  
  dat = list(name1=q$name, name2=t$name, 
    seg1=qPanel, seg2=tPanel, xaxis1=xaxis1, xaxis2=xaxis2, 
    comp=comp, 
    gap1=gap1, gap2=gap2, mapp1=mapp1, mapp2=mapp2,
    gene1=gene1, gene2=gene2,
    max_len=max_len)
}



data_preprocess_old <- function(tas, dat1, dat2) {
  dfm = assign_block_mapping(tas, 10000, 10000)$df
  dfm = dfm[order(dfm$hId, dfm$hBeg, dfm$hEnd), !colnames(dfm) %in% c('block')]

  list1 = assign_block(dfm[,c('qId', 'qBeg', 'qEnd', 'strand')])
  list2 = assign_block(dfm[,c('hId', 'hBeg', 'hEnd')])

  dfb1.1 = merge(list1$dfb, dat1$seqlen, by='id')
  dfb1.2 = prepare_coord(dfb1.1)
  dfb1 = dfb1.2

  dfb2.1 = merge(list2$dfb, dat2$seqlen, by='id')
  dfb2.2 = prepare_coord(dfb2.1)
  dfb2 = dfb2.2

  max_len = max(dfb1$end.a+dfb1$beg.a[1]-1, dfb2$end.a+dfb2$beg.a[1]-1)
  tick_itv = diff( pretty(c(1,max(dfb1$len, dfb2$len)))[1:2] )
  xaxis1 = get_ticks(dfb1, tick_itv)
  xaxis2 = get_ticks(dfb2, tick_itv)

  dfm.a = filter_assign_block(dfm[,c('qId','qBeg','qEnd','strand')], dfb1)
  dfm.b = filter_assign_block(dfm[,c('hId','hBeg','hEnd')], dfb2)
  stopifnot(rownames(dfm.a) == rownames(dfm), rownames(dfm.b) == rownames(dfm))
  dfm.2 = cbind(dfm, blk1=dfm.a$block, blk1.beg.a=dfm.a$beg.a, blk1.end.a=dfm.a$end.a, blk1.strand.a=dfm.a$strand.a, blk2=dfm.b$block, blk2.beg.a=dfm.b$beg.a, blk2.end.a=dfm.b$end.a, blk2.strand.a=dfm.b$strand.a, stringsAsFactors=F)
  comp = dfm.2
  
  name1 = dat1.name; name2 = dat2.name
  
  if(empty(dat1$gap)) {gap1=NULL} else {gap1=filter_assign_block(dat1$gap, dfb1)}
  if(empty(dat2$gap)) {gap2=NULL} else {gap2=filter_assign_block(dat2$gap, dfb2)}

  if(empty(dat1$gene)) {gene1=NULL} else {gene1=filter_assign_block(dat1$gene, dfb1)}
  if(empty(dat2$gene)) {gene2=NULL} else {gene2=filter_assign_block(dat2$gene, dfb2)}
  
  if(empty(dat1$te)) {te1=NULL} else {te1=filter_assign_block(dat1$te, dfb1)}
  if(empty(dat2$te)) {te2=NULL} else {te2=filter_assign_block(dat2$te, dfb2)}

  if(empty(dat1$nbs)) {nbs1=NULL} else {nbs1=filter_assign_block(dat1$nbs, dfb1)}
  if(empty(dat2$nbs)) {nbs2=NULL} else {nbs2=filter_assign_block(dat2$nbs, dfb2)}

  if(empty(dat1$crp)) {crp1=NULL} else {crp1=filter_assign_block(dat1$crp, dfb1)}
  if(empty(dat2$crp)) {crp2=NULL} else {crp2=filter_assign_block(dat2$crp, dfb2)}
  
  if(empty(dat1$mapp)) {mapp1=NULL} else {mapp1=filter_assign_block_wig(dat1$mapp, dfb1)}
  if(empty(dat2$mapp)) {mapp2=NULL} else {mapp2=filter_assign_block_wig(dat2$mapp, dfb2)}
  
  list(name1=name1, name2=name2, seg1=dfb1, seg2=dfb2, xaxis1=xaxis1, xaxis2=xaxis2, comp=comp, 
    gap1=gap1, gap2=gap2, mapp1=mapp1, mapp2=mapp2,
    gene1=gene1, gene2=gene2, te1=te1, te2=te2, nbs1=nbs1, nbs2=nbs2, crp1=crp1, crp2=crp2,
    max_len=max_len)
}

## plotting functions
plot_segment <- function(df, y=unit(0.5,'npc'), col.p='red', col.n='blue', text.above=F, text.rot=0, vp=NULL) {
  colnames(df) = c('id','beg','end','strand')
  idxs.n = which(df$strand == "-")
  
  line.cols = rep(col.p, nrow(df))
  line.cols[idxs.n] = col.n
  grid.segments(
    x0 = unit(df$beg, 'native'), x1 = unit(df$end, 'native'),
    y0 = y, y1 = y, 
    gp = gpar(col=line.cols), vp = vp)
  
  for (i in 1:nrow(df)) {
    beg = df[i,2]; end = df[i,3]; strand = df[i,4]
    if(strand == '-') {
      x1 = unit(beg, 'native') + unit(5,'points')
      x2 = unit(beg, 'native')
    } else {
      x1 = unit(end, 'native') - unit(5,'points')
      x2 = unit(end, 'native')
    }
    if(i == 1) {
      arrows.x1 = x1
      arrows.x2 = x2
    } else {
      arrows.x1 = unit.c(arrows.x1, x1)
      arrows.x2 = unit.c(arrows.x2, x2)
    }
  }
  grid.segments(
    x0 = arrows.x1, x1 = arrows.x2,
    y0 = y - unit(3, 'points'), y1 = y, 
    gp = gpar(col=line.cols), vp = vp)
  grid.segments(
    x0 = arrows.x1, x1 = arrows.x2,
    y0 = y + unit(3, 'points'), y1 = y, 
    gp = gpar(col=line.cols), vp = vp)
  
  if(text.above) {
    text.y = y + unit(5, 'points')
  } else {
    text.y = y - unit(5, 'points')
  }
  text.offset = ifelse(text.above, unit(-10, "points"), unit(10, 'points'))
  grid.text( label = df$id, x = unit(df$beg, 'native'), 
    y = text.y, just = c("left","center"), 
    rot = text.rot, gp = gpar(cex=0.8, fontfamily="Helvetica"), vp = vp)
}
plot_xaxis <- function(xaxis, y=unit(0,'npc'), tick.above=F, vp=NULL) {
  dfl = xaxis$line
  dft = xaxis$tick
  grid.segments(
    x0 = unit(dfl$beg.a, 'native'),
    x1 = unit(dfl$end.a, 'native'),
    y0 = y, y1 = y, vp = vp)
  
  if( tick.above ) {
    tick.y = y + unit(3, 'points')
    text.y = y + unit(5, 'points')
    text.just = c('center', 'bottom')
  } else {
    tick.y = y - unit(3, 'points')
    text.y = y - unit(5, 'points')
    text.just = c('center', 'top')
  }
  grid.segments(
    x0 = unit(dft$pos.a, 'native'),
    x1 = unit(dft$pos.a, 'native'),
    y0 = y, y1 = tick.y, 
    vp = vp)
  grid.text( 
    label = dft$pos/1000,
    x = unit(dft$pos.a, 'native'), 
    y = text.y,  just = text.just, 
    gp = gpar(cex=0.7), vp = vp)
}
plot_feature_ds <- function(df, y=unit(0.5,'npc'), height=unit(5,'points'), fill='slategray3', fill.p=fill, fill.n=fill, text.show=F, text.offset=unit(15,'points'), text.above=F, text.rot=0, vp=NULL) {
  colnames(df) = c('id','beg','end','strand')
  dfp = df[df$strand == "+",]
  dfn = df[df$strand == "-",]
  if( !empty(dfp) ) {
    grid.rect( 
    x = unit(dfp$beg, 'native'), y = y + unit(2, 'points'),
    width = unit(dfp$end-dfp$beg, 'native'), 
    height = height, just = c('left', 'bottom'),
    gp = gpar(lwd=0, fill=fill.p, alpha=0.9), vp = vp)
  }
  if( !empty(dfn) ) {
    grid.rect( 
    x = unit(dfn$beg, 'native'), y = y - unit(2, 'points'),
    width = unit(dfn$end-dfn$beg, 'native'), 
    height = height, just = c('left', 'top'),
    gp = gpar(lwd=0, fill=fill.n, alpha=0.9), vp = vp)
  }

  if( text.show ) {
    text.y = ifelse(text.above, y+text.offset, y-text.offset)
    grid.text( df$id, 
    x = unit(df$beg, 'native'), y = text.y, just = c("left","center"), 
    rot = text.rot, gp = gpar(cex=0.8, fontfamily="Helvetica-Narrow"), vp = vp)
  }
}
plot_feature <- function(df, y=unit(0.5,'npc'), height=unit(5,'points'), fill='grey', vp=NULL) {
  colnames(df) = c('id','beg','end')
  grid.rect( 
    x = unit(df$beg, 'native'), 
    y = y,
    width = unit(df$end-df$beg, 'native'), height=height,
    just=c('left','center'),
    gp = gpar(lwd=0, fill=fill, alpha=0.9), vp = vp)
}
plot_comparison <- function(comp, y1=unit(0.95, 'npc'), y2=unit(0.05, 'npc'), fill.p='skyblue1', fill.n='tomato', alpha=0.1, vp=NULL) {
  comp.xs = c()
  comp.fills = c()
  comp.ids = c()
  for (i in 1:nrow(comp) ) {
    if( comp$qPanel.strand[i] == comp$tPanel.strand[i] ) {
      comp.x = c(comp$qPanel.beg[i], comp$qPanel.end[i], comp$tPanel.end[i], comp$tPanel.beg[i])
      comp.fill = fill.p
      comp.id = rep(i, 4)
    } else {
      comp.x = c(comp$qPanel.beg[i], comp$qPanel.end[i], comp$tPanel.beg[i], comp$tPanel.end[i])
      comp.fill = fill.n
      comp.id = rep(i, 4)
    }
    comp.xs = c(comp.xs, comp.x)
    comp.fills = c(comp.fills, comp.fill)
    comp.ids = c(comp.ids, comp.id)
  }
  
  grid.polygon(
    x = unit(comp.xs, 'native'),
    y = rep(c(y1,y1,y2,y2), nrow(comp)),
    id = comp.ids,
    gp = gpar(fill=comp.fills, alpha=alpha, lwd=0),
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
plot_title <- function(main="", subtitle="", vp=NULL) {
  grid.text(main, 
    x = unit(0.5,'npc'), y = unit(1,'npc'), 
    gp = gpar(cex=1.5, fontface='bold', fontfamily='Helvetica'),
    vp=vp)
  grid.text(subtitle, 
    x = unit(0.5,'npc'), y = unit(1,'npc') - unit(2,'lines'), 
    gp = gpar(cex=1, fontface='bold', fontfamily='Helvetica'),
    vp=vp)
}
plot_legend <- function(fill, x=unit(0.5,'npc'), y=unit(0.5,'npc'), height=unit(5,'points'), width=unit(30,'points'), vp=NULL) {
  n = length(fill)
  fill.labels = names(fill)
  grid.rect( 
    x = rep(x, n), 
    y = y - unit(seq(0, by=10, length.out=n), 'points'),
    width = width, height=height,
    just=c('left','center'),
    gp = gpar(lwd=0, fill=fill, alpha=0.9), vp = vp)
  grid.text( fill.labels, 
    x = x + width + unit(10, 'points'), 
    y = y - unit(seq(0, by=10, length.out=n), 'points'), 
    just = c("left","center"), 
    gp = gpar(cex=0.9, fontfamily="mono"), vp = vp)
}
plot_scale <- function(max_len, x=unit(0.5,'npc'), y=unit(0.5,'npc'), vp=NULL) {
  len = diff( pretty(1:max_len, 20)[1:2] )
  name = sprintf("%.0fkb", len/1000)
  grid.segments( 
    x0 = x, x1 = x + unit(len, 'native'),
    y0 = y, y1 = y,
    vp = vp)
  grid.segments( 
    x0 = unit.c(x, x + unit(len, 'native')),
    x1 = unit.c(x, x + unit(len, 'native')),
    y0 = rep(y, 2), 
    y1 = rep(y + unit(3, 'points'), 2), 
    vp = vp)
  grid.text( name,
    x = x + unit(len/2, 'native'), 
    y = y + unit(5, 'points'), 
    just = c("center","bottom"), 
    gp = gpar(cex=0.9, fontfamily="Helvetica"), vp = vp)
}

comp.plot <- function(fn, dat, width=2000, height=1000, subtitle="") {
  main = sprintf("(Top) %s : %s (Bottom)", toupper(dat$name1), toupper(dat$name2))

  fill = c('assembly gap'='black', 'gene'='tan', 'TE'='slategray3', 'nbs-lrr'='forestgreen', 'crp'='dodgerblue', 'mappability'='maroon')

  seg1=dat$seg1; seg2=dat$seg2; comp=dat$comp
  xaxis1=dat$xaxis1; xaxis2=dat$xaxis2
  gap1=dat$gap1; gap2=dat$gap2
  mapp1=dat$mapp1; mapp2=dat$mapp2
  gene1=dat$gene1; gene2=dat$gene2
  te1=dat$te1; te2=dat$te2
  nbs1=dat$nbs1; nbs2=dat$nbs2
  crp1=dat$crp1; crp2=dat$crp2
  max_len=dat$max_len
  
  CairoPNG(filename = fn, width = width, height = height)
  grid.newpage()
  vpf = viewport(x=unit(0.5,"npc"), y=unit(0.5,"npc"), 
    width=unit(1,"npc")-unit(4,"lines"), height = unit(1,"npc")-unit(6,"lines"), 
    just=c("center","center"), name='all')
  pushViewport(vpf)

  vptt <- viewport(x=0.5, y=0.95, width=1, height=0.1, xscale=c(1,max_len), name='title_top')
  pushViewport(vptt); upViewport()
  vplt <- viewport(x=0.5, y=0.85, width=1, height=0.1, xscale=c(1,max_len), name='legend_top')
  pushViewport(vplt); upViewport()
  vpat <- viewport(x=0.5, y=0.75, width=1, height=0.1, xscale=c(1,max_len), name='annt_top')
  pushViewport(vpat); upViewport()
  vpxt <- viewport(x=0.5, y=0.65, width=1, height=0.1, xscale=c(1,max_len), name='axis_top')
  pushViewport(vpxt); upViewport()
  vpm <- viewport(x=0.5, y=0.5, width=1, height=0.2, xscale=c(1,max_len), name='mapping')
  pushViewport(vpm); upViewport()
  vpxb <- viewport(x=0.5, y=0.35, width=1, height=0.1, xscale=c(1,max_len), name='axis_bot')
  pushViewport(vpxb); upViewport()
  vpd1 <- viewport(x=0.5, y=0.3, width=1, height=unit(30,'points'), just='top', xscale=c(1,max_len), name='data1')
  pushViewport(vpd1); upViewport()

  # top 
#  grid.rect(gp=gpar(fill='grey', alpha=0.2, lwd=0), vp = vpxt)
  plot_segment(seg1[,c('id','beg.a','end.a','strand')], y=unit(10,'points'), text.above=T, text.rot=30, vp=vpxt)
  plot_xaxis(xaxis1, y=unit(0,'npc'), tick.above=F, vp=vpxt)
  
  # bottom
#  grid.rect(gp=gpar(fill='grey', alpha=0.2, lwd=0), vp = vpxb)
  plot_segment(seg2[,c('id','beg.a','end.a','strand')], y=unit(1,'npc')-unit(10,'points'), text.above=F, text.rot=-30, vp=vpxb)
  plot_xaxis(xaxis2, y=unit(1,'npc'), tick.above=T, vp=vpxb)

  if( !empty(gene2) ) {
    plot_feature_ds(gene2[,c('note','beg.a','end.a','strand.a')], y=unit(1,'npc')-unit(10,'points'), height=unit(5,'points'), fill=fill['gene'], text.show=F, text.offset=unit(30,'points'), text.above=F, text.rot=-40, vp=vpxb)
  }
  if( !empty(te2) ) {
    plot_feature_ds(te2[,c('note','beg.a','end.a','strand.a')], y=unit(1,'npc')-unit(10,'points'), height=unit(5,'points'), fill=fill['TE'], vp=vpxb)
  }
  if( !empty(crp2) ) {
    plot_feature_ds(crp2[,c('note','beg.a','end.a','strand.a')], y=unit(1,'npc')-unit(10,'points'), height=unit(5,'points'), fill=fill['crp'], vp=vpxb)
  }
  if( !empty(nbs2) ) {
    plot_feature_ds(nbs2[,c('note','beg.a','end.a','strand.a')], y=unit(1,'npc')-unit(10,'points'), height=unit(5,'points'), fill=fill['nbs-lrr'], vp=vpxb)
  }

  # middle
  plot_comparison(comp, y1=unit(0.95, 'npc'), y2=unit(0.05, 'npc'), alpha=0.5, vp=vpm)
  if( !empty(gap1) ) {
    plot_feature(gap1[,c('id','beg.a','end.a')], y=unit(0.95,'npc')+unit(3,'points'), height=unit(5,'points'), fill=fill['assembly gap'], vp=vpm)
  }
  if( !empty(gap2) ) {
    plot_feature(gap2[,c('id','beg.a','end.a')], y=unit(0.05,'npc')-unit(3,'points'), height=unit(5,'points'), fill=fill['assembly gap'], vp=vpm)
  }
  if( !empty(mapp2) ) {
    plot_hist(mapp2[,c('beg.a','end.a','score')], fill=fill['mappability'], vp=vpd1)
  }
  
  # misc
  plot_title(main=main, subtitle=subtitle, vptt)
  plot_legend(fill, x=unit(0.05,'npc'), y=unit(0.5,'npc'), vp=vplt)
  plot_scale(max_len, x=unit(0.9,'npc'), vp=vplt)

  dev.off()
}