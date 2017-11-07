source("br.fun.R")

dirw = file.path(Sys.getenv("misc2"), "briggs", "51.camoco")
#dirw = "/home/springer/zhoux379/scratch/briggs/51.camoco"

fi = file.path(dirw, "../36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

### for how to run camoco on a workstation - see Note

### obtain camoco results from scratch
gts = c("B73", "Mo17", "B73xMo17")
for (gt in gts) {
	tiw = spread(ti[ti$Genotype == gt, -c(3,4)], Tissue, fpkm)
	expr = t(as.matrix(tiw[,-1]))
	colnames(expr) = tiw[,1]
	gids = tiw[,1]
	ng = length(gids)

	pcc.matrix = cor(asinh(expr), method = 'pearson')
	pcc = pcc.matrix[lower.tri(pcc.matrix)]

	ii = which(lower.tri(pcc.matrix))
	colidx = as.integer((ii - 1) / nrow(pcc.matrix)) + 1
	rowidx = ii - (colidx - 1) * nrow(pcc.matrix)

	pcc[pcc == 1] = 0.999999
	pcc[pcc == -1] = -0.999999
	#pcc2 = log((1+pcc) / (1-pcc)) / 2
	pcc2 = atanh(pcc)
	coexv = (pcc2 - mean(pcc2)) / sd(pcc2)
	
	ord = order(-coexv)
	idxs = ord[1:(1*1000000)]
	dw = data.frame(g1 = gids[rowidx[idxs]], g2 = gids[colidx[idxs]], coex = coexv[idxs])
	
	fo1 = sprintf("%s/01.edges.%s.tsv", dirw, gt)
	write.table(dw, fo1, sep = "\t", row.names = F, col.names = F, quote = F)
	
	fo2 = sprintf("%s/02.%s.mcl", dirw, gt)
	cmd = sprintf("mcl %s --abc -scheme 7 -I %g -o %s", fo1, 2, fo2)
	system(cmd)
	
	fo3 = sprintf("%s/03.modules.%s.tsv", dirw, gt)
	cmd = sprintf("mcl2tsv.py %s %s", fo2, fo3)
	system(cmd)
	
	fm = sprintf("%s/03.modules.%s.tsv", dirw, gt)
	tm = read.table(fm, header = T, as.is = T, sep = "\t")
	cat(sprintf("%s: %d modules with %d genes\n", gt, length(unique(tm$V1)), nrow(tm)))
}


### look at B&M expression correlation v.s. connectivity correlation => re-wiring
gts = c("B73", "Mo17", "B73xMo17")

tissue = "root_0DAP"
#tiw = spread(ti[ti$Tissue == tissue, -c(2,4)], Genotype, fpkm)
ti2 = ti
ti2$fpkm = asinh(ti2$fpkm)
tiw = spread(ti2[,-4], Genotype, fpkm)
cor.test(tiw$B73, tiw$Mo17)
dft = ddply(tiw, .(Tissue), myfunc <- function(x) cor.test(x[,'B73'], x[,'Mo17'])$estimate)
ptitles = sprintf("%s: %.02f", dft$Tissue, dft$cor)
names(ptitles) = dft$Tissue

p1 = ggplot(tiw, aes(x = B73, y = Mo17)) +
  geom_point(size = 0.2) +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = 'asinh(B73 FPKM)') +
  scale_y_continuous(name = 'asinh(Mo17 FPKM)') +
  facet_wrap(~Tissue, ncol = 5, labeller = as_labeller(ptitles)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  #theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0,0), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = sprintf("%s/51.expr.corr.pdf", diro)
ggsave(p1, filename = fp, width = 10, height = 8)

#
fc1 = sprintf("%s/11.coex.%s.rda", diro, gts[1])
x = load(fc1)
coex1 = coex
de1 = apply(coex1, 1, myfunc <- function(x) sum(x >= 3))

fc2 = sprintf("%s/11.coex.%s.rda", diro, gts[2])
x = load(fc2)
coex2 = coex
de2 = apply(coex2, 1, myfunc <- function(x) sum(x >= 3))

cor.test(de1, de2)

ptitle = sprintf("PCC = %.03f", cor.test(de1, de2)$estimate)
tt = data.frame(b = de1, m = de2)
p2 = ggplot(tt, aes(x = b, y = m)) +
  geom_point(size = 0.5) +
  geom_smooth(method="lm") +
  scale_x_continuous(name = 'B73 Gene Degree') +
  scale_y_continuous(name = 'Mo17 Gene Degree') +
  theme_bw() +
  ggtitle(ptitle) + 
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.3,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0)) + 
  theme(plot.title = element_text(hjust = 0.5))
fp = sprintf("%s/52.degree.corr.pdf", diro)
ggsave(p2, filename = fp, width = 6, height = 6)


### compare camoco and wgcna
fc1 = sprintf("%s/11.coex.%s.rda", diro, gt)
x = load(fc1)
coex1 = coex

fc2 = file.path(dirw, "52.wgcna/x.RData")
x = load(fc2)
coex2 = TOM
tom = coex2[low.tri(coex2)]
pdf(file.path(dirw, "52.wgcna/coex.pdf"), width = 6, height = 6)
hist(tom, breaks = seq(0, 1, by = 0.025))
dev.off()

ec = sapply(1:nrow(coex1), myfunc <- function(i) cor(coex1[i, -i], coex2[i, -i]))
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

tp = data.frame(gid = ti$gid, ec = ec, novlp = novlps, stringsAsFactors = F)
tp = merge(tp, tg2, by.x = 'gid', by.y = 'par')
tp = tp[order(tp$ec),]
save(tp, file = file.path(dirw, "22.rda"))

p1 = ggplot(tp) +
  geom_point(aes(x = ec, y = novlp), size = 0.5) +
  scale_x_continuous(name = 'Expression Conservation (EC) Score') +
  scale_y_continuous(name = '# Common Neighbors btw. Networks', limits = c(0, 1000)) + 
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, color = "blue")) +
  theme(axis.text.y = element_text(size = 8, color = "blue", angle = 90, hjust  = 0.5))
fp = file.path(dirw, "51.camoco/21.ec.pdf")
ggsave(p1, filename = fp, width = 6, height = 6)

x = load(file.path(dirw, "22.rda"))
require(MASS)
p2 = ggplot(tp, aes(x = ec, y = novlp)) +
  stat_density2d(aes(fill = ..density..), contour = F, geom="tile") +
  scale_x_continuous(name = 'Expression Conservation (EC) Score', limits = c(-1, 1), expand = c(0, 0)) +
  scale_y_continuous(name = '# Common Neighbors btw. Networks', limits = c(0, 1000), expand = c(0, 0)) + 
  scale_fill_gradient(low = "white", high = "royalblue") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.85), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA), legend.key.size = unit(1, 'lines'), legend.title = element_text(size = 8, angle = 0), legend.text = element_text(size = 8, angle = 0)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, color = "blue")) +
  theme(axis.text.y = element_text(size = 8, color = "blue", angle = 90, hjust  = 0.5))
fp = file.path(dirw, "51.camoco/22.ec.density.pdf")
ggsave(p2, filename = fp, width = 6, height = 6)
