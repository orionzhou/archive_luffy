require(plyr)
require(dplyr)
require(ape)
require(ggplot2)
require(ggsignif)
require(WGCNA)
require(tidyr)
require(RColorBrewer)
require(igraph)
require(pheatmap)
source("circlePlot.R")

options(stringsAsFactors = FALSE)

dirw = '/home/springer/zhoux379/scratch/briggs2'

fl = file.path(dirw, '00.1.read.correct.tsv')
tl = read.table(fl, header = T, sep = "\t", as.is = T)[,1:5]
tissues = unique(tl$Tissue)

gts = c("B73", "Mo17", "B73xMo17")

#allowWGCNAThreads()
enableWGCNAThreads()

## make expression heatmap
fpkm_heatmap <- function(gids, ti, tissues, gts, fo) {
	ti1 = ti[ti$gid %in% gids & ti$Genotype %in% gts,-4]
	ti1$Tissue = factor(ti1$Tissue, levels = tissues)
	ti1$Genotype = factor(ti1$Genotype, levels = gts)
	ti2 = ti1[order(ti1$Genotype, ti1$Tissue),]
	ti2 = cbind(ti2, cond = sprintf("%s.%s", ti2$Genotype, ti2$Tissue))
	ti3 = ti2[,c(1,4,5)]
	ti3$cond = factor(ti3$cond, levels = unique(ti3$cond))
	ti4 = spread(ti3, cond, fpkm)
	td = ti4
	rownames(td) = td$gid
	td = td[,-1]
	td = asinh(td)
	
	ta = unique(ti2[,c(2:3,5)])
	rownames(ta) = ta$cond
	ta = ta[,-3]
	
	drows1 <- "correlation"
	dcols1 <- "correlation"
	col.pal <- brewer.pal(9, "Blues")
	col.geno = brewer.pal(8, "Paired")[6:4]
	col.tissue = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))
	names(col.tissue) = unique(ta$Tissue)
	names(col.geno) = unique(ta$Genotype)
	ann_colors = list(
		Genotype = col.geno,
		Tissue = col.tissue
	)
	
	hm.parameters <- list(td, 
		color = col.pal,
		cellwidth = 5, cellheight = 5, scale = "none",
		treeheight_row = 150,
		kmeans_k = NA,
		show_rownames = T, show_colnames = F,
		main = "Heatmap of asinh(FPKM)",
		clustering_method = "complete",
		cluster_rows = T, cluster_cols = F,
		clustering_distance_rows = drows1, 
		clustering_distance_cols = dcols1,
		annotation_col = ta,
		annotation_colors = ann_colors,
		gaps_col = c(17,34),
		fontsize_row = 6
	)
	do.call("pheatmap", c(hm.parameters, filename=fo))
}

## single-module circle plot
make_circleplot <- function() {
	pathlabel1 = sprintf("%d: Zsummary = %.02f; medianRank = %d", pid,
		tz$Mo17[tz$modName == as.character(pid)], 
		tm$Mo17[tm$modName == as.character(pid)])
	pathlabel2 = sprintf("%d: Zsummary = %.02f; medianRank = %d", pid,
		tz$B73xMo17[tz$modName == as.character(pid)], 
		tm$B73xMo17[tm$modName == as.character(pid)])

	idxGenes = which(mids == pid)
	nPathGenes = length(idxGenes)
	pathwayAdjs = list()
	KMEpathway = matrix(0, nPathGenes, nSets)
	for (set in 1:nSets)
	{
	  bc = bicor(multiExpr[[set]]$data[, idxGenes], use = "p")
	  pathwayAdjs[[set]] = abs(bc)^4 * sign(bc)
	  KMEpathway[, set] = bicor(multiExpr[[set]]$data[, idxGenes], multiMEs[[set]]$data[, mid], use = "p")
	}

	conn = matrix(0, nPathGenes, nSets)
	for (set in 1:nSets)
	  conn[, set] = apply(abs(pathwayAdjs[[set]]), 2, sum)-1
	weights = c(3,1,5,1, 3,1,5,1);
	wMat = matrix(weights, nPathGenes, nSets, byrow = TRUE)
	wconn = apply(conn * wMat, 1, sum)
	order = order(-wconn)
	# use the gene names as lables
	labels = colnames(pathwayAdjs[[ref]])
	labels = substr(labels, 7, 14)

	fo = sprintf("%s/38.circleplot/%d.pdf", dirw, pid)
	pdf(file = fo, wi=16, h=4)
	par(mfrow =c(1,4));
	par(mar = c(0.3, 0.2, 1.5, 0.2))
	for (set in 1:nSets)
	{
	  circlePlot(pathwayAdjs[[set]], labels, order, main = setLabels[set],
				variable.cex.labels = TRUE,
				radii = c(0.56, 0.62), center = c(0.1, 0.04),
				min.cex.labels = 1.2, max.cex.labels = 1.4, cex.main = 1.4)
	  if(set == 1)
		text(0,-1, pathlabel)
	}

	par(mar = c(3.3, 3.3, 4, 0.5));
	par(mgp = c(1.9, 0.6, 0))
	verboseScatterplot(KMEpathway[, ref], KMEpathway[, test],
					   xlab = sprintf("KME in %s", setLabels[ref]),
					   ylab = sprintf("KME in %s", setLabels[test]),
					   main = sprintf("%s. KME in %s vs. %s", 
						   LETTERS[1], setLabels[test], setLabels[ref]),
					   cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.2, abline = TRUE
	)
	dev.off()
}