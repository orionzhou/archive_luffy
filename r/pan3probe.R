require(GenomicRanges)
require(plyr)
require(rtracklayer)
require(parallel)

dir = '/home/youngn/zhoup/Data/misc3/pan3probe'

fg = '/home/youngn/zhoup/Data/genome/Mtruncatula_4.0/40_gene.gtb'
tg = read.table(fg, header=T, sep="\t", as.is=T, quote="")[,1:6]

gg = GRanges(seqnames=Rle(tg$chr), ranges=IRanges(tg$beg, end=tg$end))
ggr = reduce(gg)

sum(width(gg))
sum(width(ggr))

to = data.frame(id=as.character(seqnames(ggr)), beg=as.numeric(start(ggr)), end=as.numeric(end(ggr)), len=as.numeric(width(ggr)))
#write.table(to[,1:3], sprintf("%s/01_win.tbl", dir), col.names=F, row.names=F, sep='\t', quote=F)


dir = '/home/youngn/zhoup/Data/genome/pan3/18_stat_k60'

if(FALSE) {

bwg = import(file.path(dir, '11_gc.bw'), which=grgr, asRangedData=T)
gc = bwg$score
chrs = seqnames(bwg)
poss = start(bwg)
rm(bwg)

bwm = import(file.path(dir, '15_mapp.bw'), which=grgr, asRangedData=T)
mapp = bwm$score
rm(bwm)

bwd = import(file.path(dir, '17_deltag.bw'), which=grgr, asRangedData=T)
deltag = bwd$score
rm(bwd)

deltagf = deltag[deltag<900]
length(deltag)
length(deltagf)
summary(deltagf)
mean(deltagf) - 2*sd(deltagf)

gc_min = 0.3
gc_max = 0.583
mapp_min = 1
mapp_max = 1
deltag_min = -7.51
deltag_max = 20

idxs = (gc>=gc_min & gc<=gc_max & mapp>=mapp_min & mapp<=mapp_max & deltag>=deltag_min & deltag<=deltag_max)
gs = GRanges(seqnames=Rle(chrs[idxs]), ranges=IRanges(poss[idxs], end=poss[idxs]))
gsr = reduce(grs)

}

fps = "/home/youngn/zhoup/Data/misc3/pan3probe/11_win.tbl"
tps = data.frame(chr=as.character(seqnames(ggr)), beg=start(ggr), end=end(ggr))
#write.table(tps, fps, col.names=F, row.names=F, sep='\t', quote=F)

tps = read.table(fps, header=F, sep='\t', as.is=T, quote="")
colnames(tps) = c('chr', 'beg', 'end')
gps = GRanges(seqnames=Rle(tps$chr), ranges=IRanges(tps$beg, end=tps$end))

probe_dist = 150

len_ovlp <- function(x, gr) {
	require(GenomicRanges)
	chr = as.character(x['chr'])
	beg = as.numeric(x['beg'])
	end = as.numeric(x['end'])
	cat(chr, beg, end, "\n", sep=" ")
	grg1 = GRanges(seqnames=Rle(chr, 1), ranges=IRanges(beg, end))
	
	#poss = start(grs[(seq(1,n)-1) * probe_dist + 1])
	#tgp$n[i] = n
	gri = intersect(gr[seqnames(gr)==chr], grg1)
	sum(width(gri))
}

lenis = apply(tg, 1, len_ovlp, gr=gps)
tgp = cbind(tg, lenp=lenis)

tgs = tg[1:10000,]
cl = makePSOCKcluster(7)
lenis = parApply(cl, tgs, 1, len_ovlp, gr=gps)

write.table(tgp, "/home/youngn/zhoup/Data/misc3/pan3probe/12_gene.tbl", col.names=T, row.names=F, sep='\t', quote=F)
