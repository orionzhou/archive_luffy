require(plyr)
require(dplyr)
require(seqinr)
require(GenomicRanges)
require(ggplot2)
source('Location.R')

dirw = file.path(Sys.getenv('misc2'), 'pbjoin')

alg = "PBBNDT"

fs = sprintf("%s/HM340.%s/raw.fix.fas.map", Sys.getenv("genome"), alg)
ts = read.table(fs, header = F, sep = "\t", as.is = T)
colnames(ts) = c("chr", "size", "chr2")

### read in res-enzyme sites
res = "BbvCI"
f02 = sprintf("%s/02.restriction/HM340.%s.%s.txt", dirw, alg, res)
t02 = read.table(f02, header = T, sep = "\t", as.is= T)
colnames(t02) = c("chr", "motif", "beg", "end", "mm", "srd", "str")
grr1 = with(t02, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

res = "BspQI"
f02 = sprintf("%s/02.restriction/HM340.%s.%s.txt", dirw, alg, res)
t02 = read.table(f02, header = T, sep = "\t", as.is= T)
colnames(t02) = c("chr", "motif", "beg", "end", "mm", "srd", "str")
grr2 = with(t02, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

##### find joins and breaks
f01 = sprintf("%s/%s.txt", dirw, alg)
t01 = read.table(f01, header = F, sep = "\t", as.is= T)
colnames(t01) = c("nchr", "ochr", "obeg", "oend", "srd", "nbeg", "nend")

### 'joins'
x = table(t01$nchr)
tj1 = t01[t01$nchr %in% names(x)[x>1],]
tj2 = ddply(tj1, .(nchr), function(ds) {
	ds = ds[order(ds$nbeg),]
	ds = cbind(ds, type = "lr")
	ds$type = as.character(ds$type)
	ds$type[1] = 'r'
	ds$type[nrow(ds)] = 'l'
#	if(nrow(ds) > 2) {ds$type[2:(nrow(ds)-1)] = 'lr'}
	ds
})

chrs = c(); begs = c(); ends = c()
len = 10000
for (i in 1:nrow(tj2)) {
	if(tj2$type[i] %in% c('r', 'lr')) {
		chrs = c(chrs, tj2$nchr[i])
		begs = c(begs, max(tj2$nbeg[i]+1, tj2$nend[i]-len+1))
		ends = c(ends, tj2$nend[i])
	}
	if(tj2$type[i] %in% c('l', 'lr')) {
		chrs = c(chrs, tj2$nchr[i])
		begs = c(begs, tj2$nbeg[i]+1)
		ends = c(ends, min(tj2$nbeg[i]+len, tj2$nend[i]))
	}
}
tj = data.frame(chr = chrs, beg = begs, end = ends, stringsAsFactors = F)
grj = with(tj, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
grj = reduce(grj)
tj = data.frame(chr = as.character(seqnames(grj)), beg = start(grj), end = end(grj), stringsAsFactors = F)

### 'breaks'
x = table(t01$ochr)
tb1 = t01[t01$ochr %in% names(x)[x>1],]
tb1 = tb1[order(tb1$ochr, tb1$obeg),]
tb2 = ddply(tb1, .(ochr), function(ds) {
	ds = ds[order(ds$obeg),]
	ds = cbind(ds, type = "lr")
	ds$type = as.character(ds$type)
	ds$type[1] = 'r'
	ds$type[nrow(ds)] = 'l'
#	if(nrow(ds) > 2) {ds$type[2:(nrow(ds)-1)] = 'lr'}
	ds
})

chrs = c(); begs = c(); ends = c()
len = 10000
for (i in 1:nrow(tb2)) {
	if(tb2$type[i] %in% c('r', 'lr')) {
		chrs = c(chrs, tb2$nchr[i])
		begs = c(begs, max(tb2$nbeg[i]+1, tb2$nend[i]-len+1))
		ends = c(ends, tb2$nend[i])
	}
	if(tb2$type[i] %in% c('l', 'lr')) {
		chrs = c(chrs, tb2$nchr[i])
		begs = c(begs, tb2$nbeg[i]+1)
		ends = c(ends, min(tb2$nbeg[i]+len, tb2$nend[i]))
	}
}
tb = data.frame(chr = chrs, beg = begs, end = ends, stringsAsFactors = F)
grb = with(tb, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
grb = reduce(grb)
tb = data.frame(chr = as.character(seqnames(grb)), beg = start(grb), end = end(grb), stringsAsFactors = F)

tp = rbind(cbind(tj, type='join'), cbind(tb, type='break'))
tp = cbind(tp, len = tp$end-tp$beg+1)
grp = with(tp, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

cnt1 = intersect_count(grp, grr1)
cnt2 = intersect_count(grp, grr2)
tp = cbind(tp, den_bbv = cnt1/tp$len*1000, den_bsp = cnt2/tp$len*1000)

#### generate background comparison by sampling genome windows
seqlens = ts$size
names(seqlens) = ts$chr
grw = unlist(tileGenome(seqlens, tilewidth = 10000))
grw = grw[width(grw) > 5000]
tw = data.frame(chr = as.character(seqnames(grw)), beg = start(grw), end = end(grw), stringsAsFactors = F)
tw = cbind(tw, len = tw$end - tw$beg + 1)

cnt1 = intersect_count(grw, grr1)
cnt2 = intersect_count(grw, grr2)
tw = cbind(tw, den_bbv = cnt1/tw$len*1000, den_bsp = cnt2/tw$len*1000)

to1 = rbind(tp[,c('type','den_bbv')], data.frame(type='random', den_bbv=tw$den_bbv, stringsAsFactors = F))
to1 = cbind(to1, enz = 'BbvCI')
colnames(to1)[2] = 'den'
to2 = rbind(tp[,c('type','den_bsp')], data.frame(type='random', den_bsp=tw$den_bsp, stringsAsFactors = F))
to2 = cbind(to2, enz = 'BspQI')
colnames(to2)[2] = 'den'
to = rbind(to1, to2)
p1 = ggplot(to) +
  geom_density(aes(den, fill = type), alpha = 0.5, adjust = 4) + 
  scale_x_continuous(name = '# restriction enzyme cut sits (per 1kb)') +
  scale_y_continuous(name = 'Density') +
  scale_fill_brewer(palette = 'Accent') +
  facet_grid(. ~ enz) +
  theme_bw() +
  ggtitle(alg) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(legend.position = 'top', legend.background = element_blank(), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(1, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_blank(), legend.text = element_text(size = 9, angle = 0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "royalblue", angle = 0, hjust = 1))

fp = sprintf("%s/31.%s.pdf", dirw, alg)
ggsave(p1, filename = fp, width = 8, height = 5)

y1 = to$den[to$type=='join' & to$enz=='BbvCI']
y2 = to$den[to$type=='break' & to$enz=='BbvCI']
y = to$den[to$type=='random' & to$enz=='BbvCI']
t.test(y1, y, alternative = 'greater')
t.test(y2, y, alternative = 'greater')
t.test(c(y1,y2), y, alternative = 'greater')

y1 = to$den[to$type=='join' & to$enz=='BspQI']
y2 = to$den[to$type=='break' & to$enz=='BspQI']
y = to$den[to$type=='random' & to$enz=='BspQI']
t.test(y1, y, alternative = 'greater')
t.test(y2, y, alternative = 'greater')
t.test(c(y1,y2), y, alternative = 'greater')
