require(GenomicRanges)
require(dplyr)
locCluster <- function(pos, wsize) {
  npos = length(pos)
  df = data.frame(id=names(pos), pos=pos, cluster=1:npos)
  df = df[order(df$pos),]
  for (i in 1:npos) {
    for (j in (i+1):npos) {
      if(j>npos) {
        next
      }
      if(df$pos[j] - df$pos[i] <= wsize) {
        df$cluster[j] = df$cluster[i]
      }
    }
  }
  clusters = unique(df$cluster)
  tmp = cbind(cluster1=clusters, cluster2=1:length(clusters))
  x = merge(df, tmp, by.x='cluster', by.y='cluster1')
  
  df2 = data.frame(id=x$id, pos=x$pos, cluster=x$cluster2, cluster_y=1)
  df2 = df2[order(df2$pos),]
  
  clusterP = ''
  for (i in 1:npos) {
    if(df2$cluster[i] == clusterP) {
      df2$cluster_y[i] = df2$cluster_y[i-1] + 1
    } else {
      clusterP = df2$cluster[i]
    }
  }

#  hist(as.matrix(table(df2$cluster)), xlab='cluster size', main=paste(wsize, 'bp', sep=''))
  df2
}

locStr2List <- function(str) {
  tmp = gregexpr("complement([[:graph:]]+)", str)[[1]]
  strand = 1
  if(tmp[1] == 1) { strand = -1 }
  tmp = gregexpr("([0-9]+\\.\\.[0-9]+)", str)[[1]]
  strs = substring(str, tmp, tmp+attr(tmp, "match.length")-1)
  tmp = strsplit(strs, "\\.\\.")
  begs = as.numeric(sapply(tmp, "[", 1))
  ends = as.numeric(sapply(tmp, "[", 2))
  list(strand=strand, begs=begs, ends=ends)
}

locStr2Df <- function(str, seqid='chrN', type='') {
  tmp = gregexpr("complement([[:graph:]]+)", str)[[1]]
  strand = 1 
  if(tmp[1] == 1) { strand = -1 }
  tmp = gregexpr("([0-9]+\\.\\.[0-9]+)", str)[[1]]
  if(tmp[1] == -1) {
    data.frame()
  } else {
    strs = substring(str, tmp, tmp+attr(tmp, "match.length")-1)
    tmp = strsplit(strs, "\\.\\.")
    begs = as.numeric(sapply(tmp, "[", 1))
    ends = as.numeric(sapply(tmp, "[", 2))
    data.frame(chr=seqid, beg=begs, end=ends, type=type)
  }
}
get_loc_gene <- function(g) {
  id = as.character(g['id'])
  chr = as.character(g['chr'])
  locC = as.character(g['locC'])
  locI = as.character(g['locI'])
  loc5 = as.character(g['loc5'])
  loc3 = as.character(g['loc3'])
  df = rbind(locStr2Df(locC,chr,"cds"), locStr2Df(locI,chr,'intron'), locStr2Df(loc5,chr,'utr5'), locStr2Df(loc3,chr,'utr3'))
  cbind(df, id=id)
}

intersect_basepair <- function(gr1, gr2) {
  t1 = data.frame(idx = 1:length(gr1), chr = seqnames(gr1), beg = start(gr1),
    end = end(gr1), stringsAsFactors = F)
  t2 = data.frame(idx = 1:length(gr2), chr = seqnames(gr2), beg = start(gr2),
    end = end(gr2), stringsAsFactors = F)
  grl = split(gr1, 1:length(gr1))
  
  ma = as.matrix(findOverlaps(gr2, grl))
  dx1 = data.frame(tidx = ma[,2], qidx = ma[,1])
  dx2 = merge(dx1, t2, by.x = 'qidx', by.y = 'idx') 

  t11 = merge(t1, dx2, by.x = 'idx', by.y = 'tidx', all.x = T)
  gp = group_by(t11, idx)
  t12 = summarise(gp, leno = sum(pmin(end.x, end.y) - pmax(beg.x, beg.y) + 1))
  t12$leno[is.na(t12$leno)] = 0
  t12$leno
}
intersect_score <- function(gr1, gr2) {
  t1 = data.frame(idx = 1:length(gr1), chr = seqnames(gr1), beg = start(gr1),
    end = end(gr1), stringsAsFactors = F)
  t2 = data.frame(idx = 1:length(gr2), chr = seqnames(gr2), beg = start(gr2),
    end = end(gr2), score = mcols(gr2)$score, stringsAsFactors = F)
  grl = split(gr1, 1:length(gr1))
  
  ma = as.matrix(findOverlaps(gr2, grl))
  dx1 = data.frame(tidx = ma[,2], qidx = ma[,1])
  dx2 = merge(dx1, t2, by.x = 'qidx', by.y = 'idx') 

  t11 = merge(t1, dx2, by.x = 'idx', by.y = 'tidx', all.x = T)
  gp = group_by(t11, idx)
  t12 = summarise(gp, score = sum(score))
  t12$score[is.na(t12$score)] = 0
  t12$score
}
intersect_count <- function(gr1, gr2) {
  t1 = data.frame(idx = 1:length(gr1), chr = seqnames(gr1), beg = start(gr1),
    end = end(gr1), stringsAsFactors = F)
  t2 = data.frame(idx = 1:length(gr2), chr = seqnames(gr2), beg = start(gr2),
    end = end(gr2), stringsAsFactors = F)
  grl = split(gr1, 1:length(gr1))
  
  ma = as.matrix(findOverlaps(gr2, grl))
  dx1 = data.frame(tidx = ma[,2], qidx = ma[,1])
  dx2 = merge(dx1, t2, by.x = 'qidx', by.y = 'idx') 

  t11 = merge(t1, dx2, by.x = 'idx', by.y = 'tidx')
  t12 = transform(t11, leno = pmin(end.x, end.y) - pmax(beg.x, beg.y) + 1)
  t13 = t12[t12$leno / (t12$end.y - t12$beg.y + 1) >= 0.5, ]
  gp = group_by(t13, idx)
  t21 = summarise(gp, cnt = length(leno))
  
  cnts = rep(0, nrow(t1))
  names(cnts) = t1$idx
  cnts[t21$idx] = t21$cnt
  as.numeric(cnts)
}
