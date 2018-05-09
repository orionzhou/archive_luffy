source("br.fun.R")

dirg = file.path(Sys.getenv("genome"), "Zmays_v4")
dirw = '/home/springer/zhoux379/data/misc2/dev40'

## merge briggs2 and grn23
dirw1 = '/home/springer/zhoux379/scratch/briggs2'
dirw2 = '/home/springer/zhoux379/data/misc2/grn23'

fm1 = file.path(dirw1, "00.1.read.correct.tsv")
tm1 = read.table(fm1, sep = "\t", header = T, as.is = T)
sids = tm1$SampleID[tm1$Genotype == "B73"]
tm1 = tm1[tm1$SampleID %in% sids,]

fi1 = file.path(dirw1, "34.fpkm.tsv")
fi1 = file.path(dirw1, "33.fpm.tsv")
ti1 = read.table(fi1, sep = "\t", header = T, as.is = T)
ti1 = ti1[,c('gid',sids)]

fm2 = file.path(dirw2, "00.0.srr.tsv")
tm2 = read.table(fm2, sep = "\t", header = T, as.is = T)

fi2 = file.path(dirw2, "34.rpkm.tsv")
fi2 = file.path(dirw2, "33.rkm.tsv")
ti2 = read.table(fi2, sep = "\t", header = T, as.is = T)

identical(ti1$gid, ti2$gid)
tm = rbind(tm1[,1:5], tm2[,1:5])
ti = cbind(ti1, ti2[,-1])

til = gather(ti, sid, fpkm, -gid)
til = merge(til, tm[,c('SampleID','Tissue')], by.x = 'sid', by.y = 'SampleID')
grp = dplyr::group_by(til, gid, Tissue)
til2 = dplyr::summarise(grp, fpkm = mean(fpkm))
tw = spread(til2, Tissue, fpkm)
identical(tw$gid, ti1$gid)
tl = gather(tw, tissue, fpkm, -gid)

# filtering
nt = apply(tw[,-1], 1, myfun <- function(x) {sum(x>0)})
max_expr = apply(tw[,-1], 1, myfun <- function(x) max(x))

sum(nt >= 10)
sum(nt >= 10 & max_expr >= 1)
idxs = which((nt >= 10 & max_expr >= 1) | (nt >= 1 & max_expr >= 3))
length(idxs)
twf = tw[idxs,]
tlf = tl[tl$gid %in% twf$gid,]

fo = file.path(dirw, '36.filtered.wide.tsv')
write.table(twf, fo, sep = "\t", row.names = F, col.names = T, quote = F)
fo = file.path(dirw, '36.filtered.long.tsv')
write.table(tlf, fo, sep = "\t", row.names = F, col.names = T, quote = F)


### PCA
fi = file.path(dirw, '36.filtered.wide.tsv')
#fi = file.path(dirw, '34.fpkm.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

e1 = ti[,-1]
for (i in 1:ncol(e1)) {
	e1[,i] = e1[,i] / sum(e1[,i]) * 1000000
}

n_noexp = apply(e1, 1, myfunc <- function(x) sum(x<1))
n_exp = apply(e1, 1, myfunc <- function(x) sum(x>=1))

e = asinh(e1[n_noexp < 10,])
pca <- prcomp(e, center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
tp = cbind.data.frame(lab = 1:nrow(x), Tissue = rownames(x), x[,1:5], stringsAsFactors = F)
cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

p1 = ggplot(tp) +
  geom_point(aes(x = PC1, y = PC2)) +
  geom_text(aes(x = PC1, y = PC2, label = Tissue), nudge_y = 0.01) +
  scale_x_continuous(name = 'PC1 (%)') +
  scale_y_continuous(name = 'PC2 (%)') +
  #scale_color_manual(name = "", values = cols) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fp = sprintf("%s/41.qc/01.sample.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)

