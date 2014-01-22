require(GenomicRanges)
require(plyr)
require(rtracklayer)

dir = '/home/youngn/zhoup/Data/misc3/pan3probe'

fg = '/home/youngn/zhoup/Data/genome/pan3/55_noTE.gtb'
tg = read.table(fg, header=T, sep="\t", as.is=T, quote="")[,c(1,3:6,17)]
tg = cbind(tg, crp=grepl('CRP', tg$cat3), nbs=grepl('NBS', tg$cat3))

gg = GRanges(seqnames=Rle(tg$chr), ranges=IRanges(tg$beg, end=tg$end))
ggr = reduce(gg)

glg = with(tg, makeGRangesListFromFeatureFragments(seqnames=chr, fragmentStarts=sprintf("%d,",beg), fragmentWidths=sprintf("%d,",end-beg+1), strand=srd))

get_genome_stat <- function(gr) {
  dir = '/home/youngn/zhoup/Data/genome/pan3/18_stat_k60'
  bwg = import(file.path(dir, '11_gc.bw'), which=gr, asRangedData=T)
  gc = bwg$score
  chrs = seqnames(bwg)
  poss = start(bwg)
  rm(bwg)

  bwm = import(file.path(dir, '15_mapp.bw'), which=gr, asRangedData=T)
  mapp = bwm$score
  rm(bwm)

  bwd = import(file.path(dir, '17_deltag.bw'), which=gr, asRangedData=T)
  deltag = bwd$score
  rm(bwd)
  data.frame(chr=chrs, pos=poss, gc=gc, mapp=mapp, deltag=deltag)
}
pick_pos <- function(df, probe_dist) {
  poss = c(-150)
  for (i in 1:nrow(df)) {
    while( poss[length(poss)] + probe_dist < df$end[i] ) {
      bound = poss[length(poss)] + probe_dist
      if(bound < df$beg[i]) {
        poss = c(poss, df$beg[i])
      } else {
        poss = c(poss, bound+1)
      }
    }
  }
  data.frame(chr=df$chr[1], pos=poss[-1])
}
probe_dist = 150


# first-pass design of best probes (mapp=1, good GC and deltaG)
dir = '/home/youngn/zhoup/Data/genome/pan3/18_stat_k60'

ds = get_genome_stat(ggr)

deltagf = ds$deltag[ds$deltag<900]
length(deltagf)
summary(deltagf)
mean(deltagf) - 2*sd(deltagf)

gc_min = 0.3
gc_max = 0.583
mapp_min = 1
mapp_max = 1
deltag_min = -7.51
deltag_max = 20

idxs = (ds$gc>=gc_min & ds$gc<=gc_max & ds$mapp>=mapp_min & ds$mapp<=mapp_max & ds$deltag>=deltag_min & ds$deltag<=deltag_max)
gs = GRanges(seqnames=Rle(ds$chr[idxs]), ranges=IRanges(ds$pos[idxs], end=ds$pos[idxs]))
gsr = reduce(gs)

tps = data.frame(chr=as.character(seqnames(gsr)), beg=start(gsr), end=end(gsr))

tpr = ddply(tps, .(chr), pick_pos, probe_dist)
#tpr1 = cbind(tpr, end = dpr$pos+60-1, id=sprintf("MtPan3_%s_%s", dpr$chr, dpr$pos))

fpr = "/home/youngn/zhoup/Data/misc3/pan3probe/15_probes.tbl"
write.table(tpr, fpr, col.names=F, row.names=F, sep='\t', quote=F)

# work on failed genes
fpr = "/home/youngn/zhoup/Data/misc3/pan3probe/15_probes.tbl"
tpr = read.table(fpr, header=F, sep='\t', as.is=T, quote="")
colnames(tpr) = c("chr", "pos")
gpr = GRanges(seqnames=Rle(tpr$chr), ranges=IRanges(tpr$beg, end=tpr$end))

midx = as.matrix(findOverlaps(gpr, glg))
didx = data.frame(gidx=midx[,2], ridx=midx[,1])
tgs1 <- ddply(didx, .(gidx), nrow)
colnames(tgs1)[2] = "nprobe1"

# second-pass design of alternative probes (mapp>=0.5)
fpr = "/home/youngn/zhoup/Data/misc3/pan3probe/15_probes.tbl"
tpr = read.table(fpr, header=F, sep='\t', as.is=T, quote="")
colnames(tpr) = c("chr", "pos")
tpr = cbind(tpr, beg=tpr$pos-probe_dist, end=tpr$pos+probe_dist)
gpr = GRanges(seqnames=Rle(tpr$chr), ranges=IRanges(tpr$beg, end=tpr$end))
ga = setdiff(ggr, gpr)

ds = get_genome_stat(ga)

gc_min = 0.3
gc_max = 0.583
mapp_min = 0.5
mapp_max = 1
deltag_min = -7.51
deltag_max = 20

idxs = (ds$gc>=gc_min & ds$gc<=gc_max & ds$mapp>=mapp_min & ds$mapp<=mapp_max & ds$deltag>=deltag_min & ds$deltag<=deltag_max)
gs = GRanges(seqnames=Rle(ds$chr[idxs]), ranges=IRanges(ds$pos[idxs], end=ds$pos[idxs]))
gsr = reduce(gs)

tps = data.frame(chr=as.character(seqnames(gsr)), beg=start(gsr), end=end(gsr))

tpr = ddply(tps, .(chr), pick_pos, probe_dist)

fpr = "/home/youngn/zhoup/Data/misc3/pan3probe/16_probes_alt.tbl"
write.table(tpr, fpr, col.names=F, row.names=F, sep='\t', quote=F)

# stats per gene
fp1 = "/home/youngn/zhoup/Data/misc3/pan3probe/15_probes.tbl"
tp1 = read.table(fp1, header=F, sep="\t", as.is=T)
gp1 = GRanges(seqnames=Rle(tp1$V1), ranges=IRanges(tp1$V2, end=tp1$V2))

midx = as.matrix(findOverlaps(gp1, glg))
didx = data.frame(gidx=midx[,2], ridx=midx[,1])
tgs1 <- ddply(didx, .(gidx), nrow)
colnames(tgs1)[2] = "nprobe1"

fp2 = "/home/youngn/zhoup/Data/misc3/pan3probe/16_probes_alt.tbl"
tp2 = read.table(fp2, header=F, sep="\t", as.is=T)
gp2 = GRanges(seqnames=Rle(tp2$V1), ranges=IRanges(tp2$V2, end=tp2$V2))

midx = as.matrix(findOverlaps(gp2, glg))
didx = data.frame(gidx=midx[,2], ridx=midx[,1])
tgs2 <- ddply(didx, .(gidx), nrow)
colnames(tgs2)[2] = "nprobe2"

tg1 = cbind(tg, gidx=1:nrow(tg))
tg11 = merge(tg1, tgs1, by="gidx", all.x=T)
tgp = merge(tg11, tgs2, by="gidx", all.x=T)
tgp[is.na(tgp$nprobe1), 'nprobe1'] = 0
tgp[is.na(tgp$nprobe2), 'nprobe2'] = 0
tgp = cbind(tgp, nprobe=tgp$nprobe1+tgp$nprobe2)

fgp = "/home/youngn/zhoup/Data/misc3/pan3probe/21_probe_per_gene.tbl"
write.table(tgp, fgp, col.names=T, row.names=F, sep='\t', quote=F)

simplestat <- function(df) {
  table( cut(df$nprobe1, breaks=c(-1,0,1,2,Inf), labels=c(0,1,2,'>=3')) ) 
  table( data.frame( nprobe=cut(df$nprobe1, breaks=c(-1,0,1,2,Inf), labels=c(0,1,2,'>=3')),  len=cut(df$end-df$beg+1, breaks=c(0,400,Inf)) ) )
}
simplestat(tgp)
simplestat(tgp[tgp$crp,])
simplestat(tgp[tgp$nbs,])
