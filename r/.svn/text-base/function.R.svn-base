readAssembly <- function(db="mt_30", opt=1) {
  f = file.path(DIR_R, "01_conf", paste("assembly_", db, ".txt", sep=""));
  a = read.table(f, header=TRUE, sep="\t");
  if(opt == 1) {
    a = a[a$chr != 'chr0', ]
    a$chr = factor(a$chr, levels=unique(a$chr));
  }
  a
}
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, colour="red", geom=geom, width=0.2, ...)
}
rowConvert <- function(x, cs, vs) {
  y = cbind(matrix(rep(as.character(x[cs]),length(vs)), nrow=length(vs), byrow=TRUE), names(x)[vs], x[vs]);
  colnames(y) = c(names(x)[cs], 'variable', 'value')
  rownames(y) = c()
  y
}
plotParCorr <- function() {
  a <- read.table(file.path(DIR_Data, "in", "PartialCorr_10kcov.txt"),header=T,sep="\t",quote="");
  a = as.matrix(a);
  cs = c(1,8);
  vs = c(2:5, 7)
  for (i in 1:nrow(a)) {
    tmp = rowConvert(a[i,], cs, vs)
    if(i == 1) {
  	b = tmp
    } else {
  	b = rbind(b, tmp)
    }
  }
  c = data.frame(b)
  c = transform(c, DistToCentroStart=as.numeric(DistToCentroStart), value=as.numeric(value))
  d = c[c$chr == 'chr1',]
  p <- ggplot(data = d)
  p + geom_area(aes(x = DistToCentroStart,y = value, fill=variable)) +
    facet_grid(variable~.) +
    opts(title='chr1')
}
strjoin = function(x, sep=" ") paste(c(x), sep=sep, collapse=sep)
lg <- function(a, b, name1, name2, f_png) {
  png(filename=f_png, width=500, height=500, units='px');
	plot(a, b, type="p", xlab=name1, ylab=name2)
	fit = lm(b~a)
	abline(fit, col="blue")
	fit.sum = summary(fit)
	ann = paste("adjusted Rsquare = ", sprintf("%.04f", fit.sum$adj.r.squared), sep="")
	text(0.8*min(a)+0.2*max(a), 0.8*min(b)+0.2*max(b), ann, col='red')
  dev.off();
}

strconcat <- function(x) { paste(x, sep=" ", collapse=" ") }
get_mt_ids <- function(opt) {
  f_id = "../conf/acc_ids.tbl"
  idt = read.table(f_id, sep="\t", header=T, stringsAsFactors=F)
  idt$id[which(idt[,opt]==1)]
}
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

  hist(as.matrix(table(df2$cluster)), xlab='cluster size', main=paste(wsize, 'bp', sep=''))
  df2
}

read_ssp <- function(fi, quiet=T) {
  line1 = scan(fi, what=integer(0), sep="\t", nlines=1, quiet=quiet)
  n_ind = line1[1]
  n_pos = line1[2]
  poss = scan(fi, what=double(0), sep="\t", nlines=1, skip=1, quiet=quiet)
  if(n_pos != length(poss)) {
    warning(paste("line 2 not", n_pos, "positions:", length(poss)))
  }
  anc = scan(fi, what=character(0), sep="\t", nlines=1, skip=2, quiet=quiet)
  if(n_pos != length(anc)) {
    warning(paste("line 3 not", n_pos, "ancestral alleles:", length(anc)))
  }
  hasOG = sum( anc %in% c('N', '?') ) < n_pos
  inds = c()
  data = matrix(NA, nrow=n_ind, ncol=n_pos)
  for (i in 1:n_ind) {
    tmp = scan(fi, what=character(), sep="\t", nlines=1, skip=i+2, quiet=quiet)
    inds = c(inds, tmp[1])
    data[i,] = tmp[-1]
  }
  rownames(data) = inds
  if(hasOG) {
    n_ind = n_ind + 1
    inds = c('anc', inds)
    data = rbind(anc, data)
  }
  list(n_ind=n_ind, n_pos=n_pos, hasOG=hasOG, inds=inds, poss=poss, data=data)
}
get_sfs_ssp = function(ssp) {
  sfs = matrix(NA, ncol=4, nrow=ssp$n_pos)
  colnames(sfs) = c("n_N", "n_anc", "n_der", "freq_der")
  rownames(sfs) = ssp$poss
  for (i in 1:ssp$n_pos) {
    if(ssp$hasOG) {
      anc = ssp$data[1,i]
      alleles = unique(ssp$data[,i])
      t = table(ssp$data[2:ssp$n_ind,i])
      if(anc == 'N') {
        warning(paste("ancestral state of position", ssp$poss[i], "is N"))
      } else if(sum(! alleles %in% c(anc, 'N')) != 1) {
        warning(paste("position", ssp$poss[i], "not having 2 alleles"))
      } else {
        der = alleles[which(!alleles %in% c(anc,'N'))]
        if(sum(alleles=='N')>0) sfs[i,'n_N'] = t['N'] else sfs[i, 'n_N'] = 0
        sfs[i, 'n_anc'] = t[anc]
        sfs[i, 'n_der'] = t[der]
      }
    } else {
      t = table(ssp$data[,i])
      alleles = names(sort(t)) 
      if(sum(alleles != 'N') != 2) {
        warning(paste("position", ssp$poss[i], "not having 2 alleles"))
      } else {
        alleles_not_N = alleles[which(alleles != 'N')]
        anc = alleles_not_N[2]
        der = alleles_not_N[1]
        if(sum(alleles=='N')>0) sfs[i,'n_N'] = t['N'] else sfs[i,'n_N'] = 0
        sfs[i, 'n_anc'] = t[anc]
        sfs[i, 'n_der'] = t[der]
      }
    }
  }
  sfs[,'freq_der'] = sfs[,'n_der'] / (sfs[,'n_anc'] + sfs[,'n_der'])
  sfs
}
get_one_sfs = function(nts, hasOG) {
  i = 0
  if(hasOG) {
    anc = nts[1]
    alleles = as.numeric(nts[-1])
    n_ind = length(nts) - 1
    if(anc == 0) {
      stop(paste("ancestral state of snp", i, "is N"))
    }
  } else {
    alleles = as.numeric(nts)
    n_ind = length(nts)
    t = sort(table(alleles[alleles!=0]))
    anc = names(t)[-1]
  }
  n_N = sum(alleles == 0)
  states = unique(alleles[alleles!=0])
  if(length(states) != 2) {
    stop(paste("snp", i, "does not have 2 alleles"))
  }
  der = states[states!=anc]
  n_anc = sum(alleles == anc)
  n_der = sum(alleles == der) 
  c(n_N, n_anc, n_der)
}
get_bam_coverage = function(f_bam, chr, beg, end) {
  system(paste("bamDepth -i ", f_bam, " -r ", chr, ":", beg, "-", end, " > /tmp/depth", sep=""))
  depth = read.table("/tmp/depth", sep="\t", header=F, as.is=T)
  colnames(depth) = c("chr", "pos", "depth")
  cov = rep(0, end - beg + 1)
  cov[depth$pos - beg + 1] = depth$depth
  cov
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
testHWE <- function(AA,Aa,aa) {
 total = AA+Aa+aa
 PA = (AA+Aa/2)/total
 Pa = 1-PA
 AAe = total * PA^2
 Aae = total * PA * Pa * 2
 aae = total * Pa * Pa
 chsq = (AA-AAe)^2/AAe + (Aa-Aae)^2/Aae + (aa-aae)^2/aae
 1-pchisq(chsq,df=1)
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


