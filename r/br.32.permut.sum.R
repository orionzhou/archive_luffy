require(plyr)
require(ape)
require(dplyr)
require(tidyr)
require(ggplot2)
require(scales)

options(stringsAsFactors = FALSE)

dirw = file.path(Sys.getenv("misc2"), "grn23", "51.camoco")
dirw = '/home/springer/zhoux379/scratch/briggs2'
diro = file.path(dirw, '51.camoco')

fi = file.path(dirw, "36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)


tde1 = data.frame(gid = unique(ti$gid))
tde2 = data.frame(gid = unique(ti$gid))
tec = data.frame(gid = unique(ti$gid))
corrs = c()

xb = seq(-5, 5, by = 0.02)
xm <- xb[-length(xb)] + 0.5 * diff(xb)
tden = data.frame(x = rep(xm, each = length(xm)), y = rep(xm, length(xm)))


for (i in 1:101) {
fo = sprintf("%s/routs/%03d.rda", diro, i)
if(!file.exists(fo)) {
	cat(sprintf("%s not there\n", fo))
	next
}
load(fo)
#save(mat.cor, de1, de2, ec, d, file = fo)

corrs = c(corrs, mat.cor)
#cat(sprintf("%d: corr %.03f\n", i, mat.cor))

td1 = data.frame(gid = gids, de1 = de1)
tde1 = merge(tde1, td1, by = 'gid', all.x = T)
colnames(tde1)[ncol(tde1)] = sprintf("permut%d", i)
td2 = data.frame(gid = gids, de2 = de2)
tde2 = merge(tde2, td2, by = 'gid', all.x = T)
colnames(tde2)[ncol(tde2)] = sprintf("permut%d", i)

tec1 = data.frame(gid = gids, ec = ec)
tec = merge(tec, tec1, by = 'gid', all.x = T)
colnames(tec)[ncol(tec)] = sprintf("permut%d", i)
#cat(sprintf("%d: dim %d %d %d\n", i, nrow(tde1), nrow(tde2), nrow(tec)))

d$x = as.numeric(levels(d$x))[d$x]
d$y = as.numeric(levels(d$y))[d$y]
d = d[d$x > 0 & d$x <= length(xm) & d$y > 0 & d$y <= length(xm),]
d$x <- xm[d$x]
d$y <- xm[d$y]
d = d[!is.na(d$x) & !is.na(d$y),]

tden = merge(tden, d, by = c('x','y'), all.x = T)
colnames(tden)[ncol(tden)] = sprintf("permut%d", i)
#cat(sprintf("%d dim %d\n", i, nrow(tden)))
}

tden[is.na(tden)] = 0
p2 = ggplot(tden, aes(x = x, y = y, fill = permut1)) +
  geom_tile() + 
  scale_fill_gradientn(name = 'Density', colours = terrain.colors(10)) +
  scale_x_continuous(name = 'Fisher Transformed Correlation, B73', expand=c(0, 0)) +
  scale_y_continuous(name = 'Fisher Transformed Correlation, Mo17', expand=c(0, 0)) +
  theme(legend.position = "none", legend.key.size = unit(1, 'lines'), legend.background = element_rect(fill = 'white', size=0)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(panel.grid.major = element_line(color = 'black')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 9, color = "black")) +
  theme(axis.text.y = element_text(size = 9, color = "black"))

fp = sprintf("%s/61.pcc.dist.pdf", diro)
ggsave(p2, filename = fp, width = 6, height = 6)


den.diff = apply(tden, 1, myfunc <- function(x) x[3] - mean(x[-c(1:3)]))
tden2 = cbind(tden[,1:2], dendiff = den.diff)

p2 = ggplot(tden2, aes(x = x, y = y, fill = dendiff)) +
  geom_tile() + 
  #scale_fill_gradientn(name = 'Density', colours = terrain.colors(10)) +
  scale_fill_gradient2(name = '', low = "blue", mid = "white", midpoint = 0, high = "red") +
  scale_x_continuous(name = 'Fisher Transformed Correlation, B73', expand=c(0, 0)) +
  scale_y_continuous(name = 'Fisher Transformed Correlation, Mo17', expand=c(0, 0)) +
  theme(legend.position = "none", legend.key.size = unit(1, 'lines'), legend.background = element_rect(fill = 'white', size=0)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(panel.grid.major = element_line(color = 'black')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 9, color = "black")) +
  theme(axis.text.y = element_text(size = 9, color = "black"))

fp = sprintf("%s/64.pcc.diff.pdf", diro)
ggsave(p2, filename = fp, width = 6, height = 6)

fp = sprintf("%s/63.mat.cor.dist.pdf", diro)
pdf(file = fp, width = 6, height = 6)
hist(corrs[-1], 40, main = NA, xlim = c(0,1), xlab = NA)
#annotate('Actual Correlation', corrs[1])
points(x = corrs[1], y = -0.2, pch = 17, col = 'red', cex = 1.3)
dev.off()

f1 = sprintf("%s/67.ec.pdf", diro)
pdf(file = f1, width = 6, height = 4)
hist(tec$permut1, 100, xlim = c(-1,1), main = NA, xlab = 'Expression Conservation')
dev.off()

f1 = sprintf("%s/67.ec.random.pdf", diro)
ecs = apply(tec[,-c(1:2)], 1, myfunc <- function(x) mean(x))
pdf(file = f1, width = 6, height = 4)
hist(ecs, 70, xlim = c(-1,1), main = NA, xlab = 'Expression Conservation')
dev.off()


ecs = tec$permut1
ecs = (ecs - mean(ecs)) / sd(ecs)
cgids = tec$gid[ecs < -3]
fo = sprintf("%s/69.gid.ec.txt", diro)
write(cgids, fo)

find_enrichment.py --obo $data/db/go/go-basic.obo --outfile 70.go.enrich.tsv 69.gid.ec.txt $genome/Zmays_v4/61.interpro/gids.txt $genome/Zmays_v4/61.interpro/11_gene.tsv

if(FALSE) {
idxs = c(1:length(ec))[order(ec, decreasing = T)][1:100]
novlps = sapply(idxs, myfunc <- function(i) {
	idxs1 = which(coex1[i,-i] >= sort(coex1[i,-i], decreasing = T)[1000])
	idxs2 = which(coex2[i,-i] >= sort(coex2[i,-i], decreasing = T)[1000])
	sum(idxs1 %in% idxs2)
})
novlps = sapply(1:length(ec), myfunc <- function(i) {
	idxs1 = which(coex1[i,-i] >= sort(coex1[i,-i], decreasing = T)[1000])
	idxs2 = which(coex2[i,-i] >= sort(coex2[i,-i], decreasing = T)[1000])
	sum(idxs1 %in% idxs2)
})
cor.test(ec, novlps)
}