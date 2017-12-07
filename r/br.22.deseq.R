require(plyr)
require(ape)
require(tidyr)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(DESeq2)
source("common.R")

dirw = file.path(Sys.getenv("misc2"), "briggs", "43.deseq")

fi = file.path(dirw, '../00.1.read.correct.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)[,1:5]

cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

### run DE-Seq
fc1 = file.path(dirw, "../32.rc.tsv")
tc1 = read.table(fc1, header = T, sep = "\t", as.is = T)
fc2 = file.path(dirw, "../32.rc.as.tsv")
tc2 = read.table(fc2, header = T, sep = "\t", as.is = T)
rownames(tc1) = tc1$gid
rownames(tc2) = tc2$gid

to = data.frame()
comps = c("B73 vs Mo17", "B73 vs B73xMo17", "Mo17 vs B73xMo17")

for (tissue in unique(ti$Tissue)) {
	#tissue = 'root_0DAP'
	for (comp in comps) {
		gts_comp = unlist(strsplit(comp, split = " vs "))
		tis = ti[ti$Tissue == tissue & ti$Genotype %in% gts_comp,]
		rownames(tis) <- tis$SampleID
		cnts1 = tc1[,tis$SampleID]
		cnts2 = tc2[,tis$SampleID]

## DE sense
dds1 = DESeqDataSetFromMatrix(countData = cnts1, colData = tis, design = ~ Genotype)
#dds1$Genotype <- factor(dds1$Genotype, levels = levels(tis$Genotype))
#dds1 = estimateSizeFactors(dds1)
deseq1 = DESeq(dds1)

#rld <- rlog(deseq1, blind=FALSE)
#fp = sprintf("%s/stats/61.pca.%s.pdf", dirw, tissue)
#pdf(fp, width=6, height=6)
#plotPCA(rld, intgroup = "Genotype", ntop = 10000)
#dev.off()

res1 = as.data.frame(results(deseq1, pAdjustMethod = "fdr"))
is.de1 = ifelse(res1$padj < .05 & res1$log2FoldChange > 1 , "up", ifelse(res1$padj < .05 & res1$log2FoldChange < -1 , "down", 0))

cnts1[rownames(cnts1) %in% rownames(res1)[is.de1 == 'up'],][1:20,]
cnts1[rownames(cnts1) %in% rownames(res1)[is.de1 == 'down'],][1:20,]

		idxs = which(is.de1 %in% c('up','down'))
		tos = data.frame(tissue = tissue, comp = comp, gid = rownames(res1)[idxs], 
			is.de = is.de1[idxs], stringsAsFactors = F)
		to = rbind(to, tos)
	}
}

fo = file.path(dirw, '01.de.tsv')
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)


### DE sharing among tissue
ff = file.path(dirw, '01.de.tsv')
tf = read.table(ff, header = T, sep = "\t", as.is = T)
tf = tf[tf$comp == 'B73 vs Mo17', -2]

grp = dplyr::group_by(tf, gid)
tfs = dplyr::summarise(grp, ntissue = length(tissue), ndirection=length(unique(is.de)))
table(tfs$ndirection)

tp = data.frame(table(tfs$ntissue))
colnames(tp) = c("nTissue", "nGene")
p1 = ggplot(tp) +
  geom_histogram(aes(x = nTissue, y = nGene), width = 0.6, stat = 'identity') +  	
  #scale_x_continuous(name = 'DOA: (F1-MP)/(HP-MP)', limits=c(-8,8)) +
  scale_y_continuous(name = 'Num DE Genes') +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fo = sprintf("%s/11.sharing.pdf", dirw)
ggsave(p1, filename = fo, width = 3, height = 3)

tk = tf
tk$is.de[tk$is.de=='up'] = 1
tk$is.de[tk$is.de=='down'] = -1
tk$is.de = as.integer(tk$is.de)
tk = spread(tk, tissue, is.de)
tk[is.na(tk)] = 0

hcl = hclust(dist(tk[,-1]), method = "ward.D")
gidsO = tk$gid[hcl$order]

tkl = gather(tk, tissue, is.de, -gid)
tkl$gid = factor(tkl$gid, levels = gidsO)
tkl$is.de = as.character(tkl$is.de)
tkl$tissue = factor(tkl$tissue, levels = unique(ti$Tissue))

p1 = ggplot(tkl) +
  geom_tile(aes(x = tissue, y = gid, fill = is.de)) + 
  #scale_x_discrete(name = '') +
  scale_y_discrete(name = 'Genes') +
  scale_fill_manual(values = c("white", "#E41A1C", "#4DAF4A"), labels = c("Not DE", "B < M", "B > M")) + 
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill='grey', size=1)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 70, hjust = 1)) +
  theme(axis.text.y = element_blank())
fo = sprintf("%s/12.sharing.heatmap.pdf", dirw)
ggsave(p1, filename = fo, width = 4.5, height = 12)



