require(plyr)
require(ape)
require(tidyr)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(GenomicRanges)

dirg = file.path(Sys.getenv("genome"), "Zmays_v4")
dirw = '/home/springer/zhoux379/scratch/briggs'
dirw = file.path(Sys.getenv("misc2"), "briggs")
diro = file.path(dirw, "41.qc")

fm = file.path(dirw, "00.1.read.correct.tsv")
tm = read.table(fm, sep = "\t", header = T, as.is = T)[,1:5]

### look at btw-replicate correlation
#fi = file.path(dirw, '33.fpm.tsv')
fi = file.path(dirw, '34.fpkm.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

e1 = ti[,-1]
stopifnot(identical(tm$SampleID, colnames(e1)))
#colnames(e1) = sprintf("%s_%s_%s_%d", tm$Sample, tm$Tissue, tm$Genotype, tm$Treatment)
dim(e1)

n_noexp = apply(e1, 1, myfunc <- function(x) sum(x<1))
e = e1[n_noexp < 50,]
dim(e)

tm2 = unique(tm[,c("Tissue", "Genotype")])
to = data.frame()
for (i in 1:nrow(tm2)) {
	tiss = tm2$Tissue[i]
	geno = tm2$Genotype[i]
	sids = tm$SampleID[tm$Tissue == tiss & tm$Genotype == geno]
	stopifnot(length(sids) == 3)
	sid_pairs = combn(sids, m = 2)
	tc = data.frame(t(sid_pairs), stringsAsFactors = F)
	colnames(tc) = c('sid1', 'sid2')
	pccs = c()
	for (j in 1:nrow(tc)) {
		pcc = cor(e[,tc$sid1[j]], e[,tc$sid2[j]], method = 'spearman')
		pccs = c(pccs, pcc)
	}
	tc = cbind(tissue = tiss, genotype = geno, tc, pcc = pccs)
	to = rbind(to, tc)
}
tp = cbind(to, x = sprintf("%s:%s", to$tissue, to$genotype))
tp$x = factor(tp$x, levels = unique(tp$x))

p1 = ggplot(tp) +
  geom_point(aes(x = pcc, y = x), shape = 4) +
  scale_x_continuous(name = 'Spearman Correlation of FPKM btw. Replicates', limits = c(0.5, 1)) +
  #scale_y_continuous() +
  #scale_color_manual(name = "", values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(legend.position = c(0.4, 0.8), legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = sprintf("%s/06.fpkm.cor.s.pdf", diro)
ggsave(p1, filename = fp, width = 6, height = 8)

### hclust
fi = file.path(dirw, '33.fpm.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

e1 = ti[,-1]
stopifnot(identical(tm$SampleID, colnames(e1)))
colnames(e1) = sprintf("%s_%s_%s_%d", tm$Sample, tm$Tissue, tm$Genotype, tm$Treatment)

n_noexp = apply(e1, 1, myfunc <- function(x) sum(x<1))
#e5 <- t(scale(t(e4), center=TRUE, scale=TRUE))
#d5 <- data.frame(d4[, c(1,2)], e5)

e = e1[n_noexp < 50,]
cor_opts = c("spearman", "pearson")
hc_opts = c("ward.D")

cor_opt = "spearman"
hc_opt = "ward.D"
for (cor_opt in cor_opts) {
for (hc_opt in hc_opts) {

plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
#e.r.dist <- cordist(e)
#e.r.hc <- hclust(e.r.dist, method='ward')
#e.r.dendro <- as.dendrogram(e.r.hc)

hc = e.c.hc
fo = sprintf("%s/01.hc.%s.%s.pdf", diro, cor_opt, hc_opt)
pdf(fo, width = 6, height = 12)
#plot(as.dendrogram(e.c.hc, hang = 0.02), cex = 0.6, ann = T, horiz = T)
plot(as.phylo(e.c.hc), cex = 0.5, label.offset = 0.02, no.margin = T)
text(0.001, 125, plot_title, adj = 0)
dev.off()

}
}

### PCA
fi = file.path(dirw, '33.fpm.tsv')
#fi = file.path(dirw, '34.fpkm.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

e1 = ti[,-1]
stopifnot(identical(tm$SampleID, colnames(e1)))
colnames(e1) = sprintf("%s_%s_%s_%d", tm$Sample, tm$Tissue, tm$Genotype, tm$Treatment)
colnames(e1) = tm$Sample

n_noexp = apply(e1, 1, myfunc <- function(x) sum(x<1))
n_exp = apply(e1, 1, myfunc <- function(x) sum(x>=1))

e = asinh(e1[n_noexp < 50,])
pca <- prcomp(e, center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
tp = cbind.data.frame(SampleID = rownames(x), x[,1:5], stringsAsFactors = F)
tp2 = merge(tp, tm[,c(1,3:5)], by = 'SampleID')
tp2$Tissue = factor(tp2$Tissue, levels = unique(tp2$Tissue))
tp2$Genotype = factor(tp2$Genotype, levels = unique(tp2$Genotype))

cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

p1 = ggplot(tp2) +
  geom_point(aes(x = PC1, y = PC2, shape = Genotype, color = Tissue)) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(name = ylab) +
  scale_color_manual(name = "", values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.3,0.1,0.1,0.1), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fp = sprintf("%s/01.sample.pca.pdf", diro)
ggsave(p1, filename = fp, width = 8, height = 7)


### LDA
require(MASS)

tissue_map = tm$Tissue; names(tissue_map) = tm$SampleID
geno_map = tm$Genotype; names(geno_map) = tm$SampleID

e1 = ti[,-1]
n_noexp = apply(e1, 1, myfunc <- function(x) sum(x<1))
idxs = which(n_noexp < 30)
length(idxs)

tl = asinh(t(ti[idxs,-1]))
colnames(tl) = ti$gid[idxs]
tl = data.frame(Tissue = tissue_map[rownames(tl)], Genotype = geno_map[rownames(tl)], tl, stringsAsFactors = F)
tl[1:5,1:5]

r <- lda(formula = Tissue ~ ., data = tl[,-2])#, prior = c(1,1,1)/3)
r$svd^2/sum(r$svd^2)

r$prior
r$counts
#r$means
#r$scaling
prop.lda = r$svd^2/sum(r$svd^2)
prop.lda

plda <- predict(object = r, newdata = tl)
tp = data.frame(Tissue = tl$Tissue, Genotype = tl$Genotype, lda = plda$x)
tp$Tissue = factor(tp$Tissue, levels = unique(tp$Tissue))
tp$Genotype = factor(tp$Genotype, levels = unique(tp$Genotype))

xlab = sprintf("LD1 (%.01f%%)", prop.lda[1]*100)
ylab = sprintf("LD2 (%.01f%%)", prop.lda[2]*100)
cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

p1 <- ggplot(tp) + 
  geom_point(aes(x = lda.LD1, y = lda.LD3, shape = Genotype, color = Tissue)) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(name = ylab) +
  scale_color_manual(name = "", values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.3,0.1,0.1,0.1), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fp = sprintf("%s/01.sample.lda.pdf", diro)
ggsave(p1, filename = fp, width = 8, height = 7)


### FPM distribution
fr = file.path(dirw, "33.fpm.tsv")
tr = read.table(fr, sep = "\t", header = T, as.is = T)
#trl = reshape(tr, direction = 'long', varying = list(2:ncol(tr)), idvar = c("gid"), timevar = "sid", v.names = 'fpm', times = colnames(tr)[2:ncol(tr)])
trl = gather(tr, sid, fpm, -gid)
trl2 = merge(trl, tm[,c(1,3:5)], by.x = 'sid', by.y = 'SampleID')
grp = dplyr::group_by(trl2, gid, Tissue, Genotype)
trl3 = dplyr::summarise(grp, fpm = mean(fpm))
trl4 = cbind(trl3, sid = sprintf("%s-%s", trl3$Tissue, trl3$Genotype))
trl4$sid = factor(trl4$sid, levels = unique(sprintf("%s-%s", tm$Tissue, tm$Genotype)))
trl4 = trl4[order(trl4$gid, trl4$sid), c(1,5,4)]

summary(trl$fpm)
p1 = ggplot(trl) +
  geom_boxplot(aes(x = sid, y = fpm), outlier.shape = NA) + #, draw_quantiles = c(0.25, 0.5, 0.75)) + 
  coord_flip() +
  scale_x_discrete(name = '', breaks = tm$SampleID, labels = sprintf("%s|%s", tm$Tissue, tm$Genotype)) +
  scale_y_continuous(name = 'FPM', limits = c(0, 30)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, color = "black", angle = 0, hjust = 1))
fp = sprintf("%s/11.fpm.pdf", diro)
ggsave(p1, filename = fp, width = 6, height = 12)

y = apply(tr[,2:ncol(tr)], 2, myfun <- function(x) {
	x1 = sort(x, decreasing = T)
	c('top[01-05]' = sum(x1[1:5]),
		'top[06-10]' = sum(x1[6:10]),
		'top[11-15]' = sum(x1[11:15]),
		'top[16-20]' = sum(x1[16:20]),
		'top[21-30]' = sum(x1[21:30]),
		'top[31-40]' = sum(x1[31:40]),
		'top[41-50]' = sum(x1[41:50])
	)
})
y = cbind(tag = rownames(y), as.data.frame(y))
yl = reshape(y, direction = 'long', varying = list(2:ncol(y)), idvar = c("tag"), timevar = "sid", v.names = 'fpm', times = colnames(y)[2:ncol(y)])

yl$tag = factor(yl$tag, levels = rev(sort(unique(yl$tag))))
p1 = ggplot(yl) +
  geom_bar(aes(x = sid, y = fpm/1000000, fill = tag), position = 'stack', stat = 'identity', width = 0.7) +
  coord_flip() +
  scale_x_discrete(name = '', breaks = tm$SampleID, labels = sprintf("%s|%s", tm$Tissue, tm$Genotype)) +
  scale_y_continuous(name = 'Proportion Reads', expand = c(0,0), limits = c(0, 1)) +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(1,1,0.1,0.1), "lines")) +
  theme(legend.position = c(0.7, 0.7), legend.direction = "vertical", legend.justification = c(0,0), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, color = "black", angle = 0, hjust = 1))
fp = sprintf("%s/11.fpm.top50.pdf", diro)
ggsave(p1, filename = fp, width = 6, height = 12)

## get top10 expressed genes w. functions
fg = file.path(dirg, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,16:18)]
gb = group_by(tg, par)
tg2 = summarise(gb, fam = names(sort(table(cat3), decreasing = T))[1])

z = ddply(trl, .(sid), myfun <- function(x) {
	x1 = x[order(x[,'fpm'], decreasing = T),][1:10,]
	x2 = cbind.data.frame(x1, crpm = as.numeric(cumsum(x1[,'fpm'])))
	x2
})
z2 = merge(z, tg2, by.x = 'gid', by.y = 'par')
z3 = merge(z2, tm[,c(1,3)], by.x = 'sid', by.y = 'SampleID')
z4 = z3[order(z3$sid, -z3$fpm),]

fo = file.path(diro, '11.fpm.top10.tsv')
write.table(z4, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### FPM correction
trw = spread(trl4, sid, fpm)
yt = ddply(trl4, .(sid), summarise, fpm = sum(fpm))
y = apply(trw[,2:ncol(trw)], 2, myfun <- function(x) {
	x1 = sort(x, decreasing = T)
	c(
		'top0' = 0,
		'top5' = sum(x1[1:5]),
		'top10' = sum(x1[1:10]),
		'top20' = sum(x1[1:20]),
		'top30' = sum(x1[1:30]),
		'top50' = sum(x1[1:50]),
		'top100' = sum(x1[1:100])
	)
})
y = cbind.data.frame(opt = rownames(y), y, stringsAsFactors = F)
#yl2 = reshape(y, direction = 'long', varying = list(2:ncol(y)), idvar = c("opt"), timevar = "sid", v.names = 'fpm', times = colnames(y)[2:ncol(y)])
yl = gather(y, sid, fpm, -opt)

### check housekeeping genes
fh = '/home/springer/zhoux379/data/genome/Zmays_v4/housekeeping/11.tsv'
thk = read.table(fh, sep = "\t", header = T, as.is = T)
hgids = thk$ngid
th = trl4[trl4$gid %in% hgids,]
th$gid = factor(th$gid, levels = thk$ngid)

opts = c("top0", "top5", "top10", "top20", "top30", "top50", "top100")
opt = opts[6]
z = yl[yl$opt == opt,]
z2 = merge(z, yt, by = 'sid')
z3 = cbind(z2[,c('sid','opt')], sf = z2$fpm.y / (z2$fpm.y - z2$fpm.x))
th2 = merge(th, z3, by = 'sid')
th3 = cbind(th2, fpmc = th2$fpm * th2$sf)

cols = brewer.pal(nrow(thk), 'Set1')
p1 = ggplot(th3) +
  geom_bar(aes(x = sid, y = fpmc, fill = gid), position = 'stack', stat = 'identity', width = 0.7) +
  coord_flip() +
  scale_x_discrete(name = '') +
  scale_y_continuous(name = 'FPM', expand = c(0,0)) +
  scale_fill_manual(values = cols, breaks = thk$ngid, labels = thk$name) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,2,0.1,0.1), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0,0), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))
fp = sprintf("%s/13.hk.%s.pdf", diro, opt)
ggsave(p1, filename = fp, width = 6, height = 8)

tb = data.frame()
for (opt in opts) {
	z = yl[yl$opt == opt,]
	z2 = merge(z, yt, by = 'sid')
	z3 = cbind(z2[,c('sid','opt')], sf = z2$fpm.y / (z2$fpm.y - z2$fpm.x))
	th2 = merge(th, z3, by = 'sid')
	th3 = cbind(th2, fpmc = th2$fpm * th2$sf)
	tb = rbind(tb, th3[,c('sid','gid','opt','fpmc')])
}
tb2 = ddply(tb, .(opt, gid), summarize, cv = (sd(fpmc)/mean(fpmc))*100)
tb2$opt = factor(tb2$opt, levels = opts)

cols = brewer.pal(length(opts), 'Paired')
p1 = ggplot(tb2) +
  geom_bar(aes(x = gid, y = cv, fill = opt), position = 'dodge', stat = 'identity', width = 0.8) +
  coord_flip() +
  scale_x_discrete(name = '', breaks = thk$ngid, labels = thk$name) +
  scale_y_continuous(name = 'C.V. of FPM across tissues') +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,2,0.1,0.1), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0,0), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))
fp = sprintf("%s/14.hk.pdf", diro)
ggsave(p1, filename = fp, width = 7, height = 7)


### FPKM distribution (obsolete)
til = reshape(ti, direction = 'long', varying = list(2:ncol(ti)), idvar = c("gid"), timevar = "sid", v.names = 'fpkm', times = colnames(ti)[2:ncol(ti)])

summary(til$rpkm)
p1 = ggplot(til) +
  geom_boxplot(aes(x = sid, y = fpkm), outlier.shape = NA) + #, draw_quantiles = c(0.25, 0.5, 0.75)) + 
  coord_flip() +
  scale_x_discrete(name = '', breaks = tm$SampleID, labels = tm$Tissue) +
  scale_y_continuous(name = 'RPKM', limits = c(0, 30)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))
fp = sprintf("%s/12.rpkm.pdf", diro)
ggsave(p1, filename = fp, width = 6, height = 9)

### Averaging replicates & compute FPKM
trl2 = merge(trl, tm[, c("SampleID", "Tissue", "Genotype")], by.x = 'sid', by.y = 'SampleID')
grp = group_by(trl2, gid, Tissue, Genotype)
trl3 = as.data.frame(summarise(grp, fpm = mean(fpm)))
#trw = reshape(trl3, direction = 'wide', timevar = c('Tissue'), idvar = c('gid'))
#colnames(trw)[2:ncol(trw)] = gsub("fpm.", "", colnames(trw)[2:ncol(trw)])

### compute FPKM
dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
f_gtb = file.path(dirg, "51.gtb")
f_tbl = file.path(dirg, "51.tbl")
tg = read.table(f_gtb, sep = "\t", header = T, as.is = T)[,1:2]
tt = read.table(f_tbl, sep = "\t", header = F, as.is = T)
colnames(tg) = c("tid", "gid")
colnames(tt) = c("chr", "beg", "end", "srd", "tid", "type", "fam")
tt2 = tt[tt$type %in% c('cds', 'utr5', 'utr3'),]
tg2 = merge(tg, tt2, by = 'tid')

gr = with(tg2, GRanges(seqnames = chr, ranges = IRanges(beg, end = end), gid = gid))
x = unlist(reduce(split(gr, elementMetadata(gr)$gid)))
tr = data.frame(gid = names(x), chr = seqnames(x), beg = start(x), end = end(x), stringsAsFactors = F)
grp = dplyr::group_by(tr, gid)
tr2 = dplyr::summarise(grp, len = sum(end - beg + 1))

to = trl3
to2 = merge(to, tr2, by = 'gid')
stopifnot(nrow(to) == nrow(to2))
to3 = cbind(to2, fpkm = to2$fpm / (to2$len / 1000))
to4 = to3[,-5]
dim(to4)
fo = file.path(dirw, '35.long.tsv')
write.table(to4, fo, sep = "\t", row.names = F, col.names = T, quote = F)


### filtering
fi = file.path(dirw, "35.long.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

grp = dplyr::group_by(ti, gid, Genotype)
ti2 = dplyr::summarise(grp, nexp = sum(fpkm>=1))
ti3 = spread(ti2, Genotype, nexp)
gids = ti3$gid[ti3$B73 >= 1 & ti3$Mo17 >= 1]
length(gids)

to = ti[ti$gid %in% gids,]
fo = file.path(dirw, '36.long.filtered.tsv')
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)



