source("br.fun.R")

dirw = file.path(Sys.getenv("misc2"), "briggs", "48.permut")
dirw = "/home/springer/zhoux379/scratch/briggs/48.permut"

fi = file.path(dirw, "../36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

gts = c("B73", "Mo17")

ti2 = ti[ti$Genotype %in% gts,]
ti3 = cbind(ti2[,-c(2:4)], treat = sprintf("%s.%s", ti2$Genotype, ti2$Tissue))


args <- commandArgs(TRUE)
i = as.numeric(args[1])


xb = seq(-5, 5, by = 0.02)
xm <- xb[-length(xb)] + 0.5 * diff(xb)

ti4 = spread(ti3, treat, fpkm)
if(i > 1) {
	treats = colnames(ti4)[-1]
	tmp = strsplit(treats, , split = "[.]")
	geno = sapply(tmp, "[", 1)
	tiss = sapply(tmp, "[", 2)
	for (j in 1:length(unique(tiss))) {
		idxs = which(tiss == tiss[j])
		geno[idxs] = sample(geno[idxs])
	}
	colnames(ti4)[-1] = sprintf("%s.%s", geno, tiss)
}
ti5 = gather(ti4, treat, fpkm, -gid)
tmp = strsplit(ti5$treat, , split = "[.]")
geno = sapply(tmp, "[", 1)
tiss = sapply(tmp, "[", 2)
ti6 = cbind(ti5[,-2], Genotype = geno, Tissue = tiss)

grp = dplyr::group_by(ti6, gid, Genotype)
ti7 = dplyr::summarise(grp, std = sd(fpkm))
ti8 = spread(ti7, Genotype, std)
gids = ti8$gid[ti8$B73 > 0 & ti8$Mo17 > 0]

for (gt in gts) {
tw = ti6[ti6$gid %in% gids,]
tw = spread(tw[tw$Genotype == gt, -3], Tissue, fpkm)

expr = t(as.matrix(tw[,-1]))
colnames(expr) = tw[,1]

pcc.matrix = cor(asinh(expr), method = 'pearson')
pcc = pcc.matrix[lower.tri(pcc.matrix)]
pcc[pcc == 1] = 0.999999
pcc[pcc == -1] = -0.999999
#pcc2 = log((1+pcc) / (1-pcc)) / 2
pcc2 = atanh(pcc)
pcc3 = (pcc2 - mean(pcc2, na.rm = T)) / sd(pcc2, na.rm = T)
head(pcc3)

ng = nrow(tw)
coex <- matrix(rep(0, ng*ng), nrow=ng)
coex[lower.tri(coex)] = pcc3
coex = t(coex)
coex[lower.tri(coex)] = pcc3

if(gt == 'B73') {
	coexb = coex
	pccb = pcc3
} else if(gt == 'Mo17') {
	coexm = coex
	pccm = pcc3
}
}

mat.cor = cor(pccb, pccm, use = 'complete.obs')

de1 = apply(coexb, 1, myfunc <- function(x) sum(x >= 3, na.rm = T))
de2 = apply(coexm, 1, myfunc <- function(x) sum(x >= 3, na.rm = T))
#cor.test(de1, de2)

ec = sapply(1:ng, myfunc <- function(i) cor(coexb[i, -i], coexm[i, -i]))

binxy <- data.frame(x=findInterval(pccb, xb), y=findInterval(pccm, xb))
d1 <- as.data.frame.table(table(binxy))
d = d1
#d$x = as.numeric(levels(d$x))[d$x]
#d$y = as.numeric(levels(d$y))[d$y]
#d$x <- xm[d$x]
#d$y <- xm[d$y]
#d = d[!is.na(d$x) & !is.na(d$y),]

fo = sprintf("%s/01.routs/%03d.rda", dirw, i)
save(gids, mat.cor, de1, de2, ec, d, file = fo)

