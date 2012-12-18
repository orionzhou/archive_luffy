dir = file.path(DIR_Misc2, "nbs/20_Mt")
f1 = file.path(dir, "01.tbl")
d1 = read.table(f1, head=TRUE, sep="\t", quote="\"")

d2 = unique(d1)
rownames(d2) = d2[,1]
d2 = d2[,-1]

lowExpFilter <- function(x, thresh) {sum(x < thresh) == length(x)}
idx = apply(d2, 1, lowExpFilter, thresh=1)

write.table(d1[idx, ], file.path(dir, "11_filtered_genes.txt"), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
d3 = d2[!idx, ]

#scale data by row for hclust
d4 <- t(scale(t(d3), center=TRUE, scale=TRUE))

#hclust
cordist <- function(x) as.dist(1-cor(t(x), method="pearson"));
d = d4

#cl = hclust(dist(d, method="euclidean"), method='ward');
cl = hclust(cordist(d), method='ward')

png(filename=file.path(dir, "23.png"), width=800, height=800, units='px')
heatmap(as.matrix(cordist(t(d))), symm=T, distfun=cordist, labCol=F, margins=c(4,12), cexRow=1.2)
dev.off()

cl.dend <- as.dendrogram(cl)
cl.group <- cutree(cl, k=10)
write(hc2Newick(cl),file=file.path(dir, '21.newick'))

#plot dendrogram and cut trees
k = 10
png(file.path(dir, "31_tree.png"), width=800, height=400, units='px')
plot(cl, hang=0.1, cex=0.1, ann=FALSE, labels=FALSE)
groups <- cutree(cl, k=k)
groups.t <- table(groups)[unique(groups[cl$order])]
m <- c(0, cumsum(groups.t))
for (i in 1:k) {
  rect(m[i] + 0.66, par("usr")[3L]+5, m[i+1]+0.33, mean(rev(cl$height)[(k-1):k]), border=2)
  text(mean(m[i:(i+1)]), mean(c(par("usr")[3L], rev(cl$height)[k])), names(m)[i+1], col='red', cex=1.5)
}
dev.off()

et = data.frame(as.table(d))
colnames(et) = c("id", "tissue", "expression")

mapping1 = data.frame(idx=cl$order, treeIdx=1:nrow(d))
mapping2 = cbind(mapping1[order(mapping1$idx),], id=rownames(d3))
et2 = merge(et, mapping2, by.x="id", by.y="id")
fig_height = 3
if( length(unique(et2$tissue)) >= 20 ) { fig_height = 7 }
p <- ggplot(et2, aes(treeIdx, tissue)) + geom_tile(aes(fill=expression)) + 
  scale_fill_gradient(low='blue', high='red') +
  scale_y_discrete(name='') +
  scale_x_continuous(name='Tree Index') +
  opts(axis.text.x = theme_text(hjust=1, size=10)) +
  opts(axis.text.y = theme_text(hjust=1, size=10))
ggsave(file.path(dir, "32_heatmap.png"), p, width=10, height=fig_height)

sum = data.frame(mapping2, group=cl.group, d)
write.table(sum, file.path(dir, "41_sum.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


#kmeans
dirW = file.path(dirO, "cl41")
cl = kmeans(d, 10)
sum = data.frame(group = cl$cluster, id=rownames(d), d)
write.table(sum, file.path(dirW, "10_sum.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
