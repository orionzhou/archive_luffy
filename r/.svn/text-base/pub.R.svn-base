dirI = file.path(DIR_Misc2, "mtgea")
dirO = file.path(DIR_R, "mtgea")
f1 = file.path(dirI, "ME1_hybridizations.txt")
meta1 = read.table(f1, head=TRUE, sep="\t", quote="\"")
f2 = file.path(dirI, "ME1_RMA_medicago.txt")
d1 = read.table(f2, sep="\t", quote="\"")
f3 = file.path(dirI, "ME1_mas_data_medicago.txt")
dm = read.table(f3, head=TRUE, sep="\t", quote="\"")
rownames(dm) = dm[,1]
dm = dm[,-1]

idx_f1 = grep("^Mtr\\.", rownames(d1))
filter1 = rep(0, nrow(d1))
filter1[c(idx_f1)] = 1

detection_call <- function(x) ! sum(x=='A') == length(x)
filter2 = apply(dm[,c(1:54)*5-1], 1, detection_call)

filter = filter1 & filter2
d2 = d1[which(filter), ]

which(duplicated(rownames(d2)))

boxplot(as.matrix(d2))

group = c(rep(1,6), rep(2,12), rep(1,36))
group1 = which(group == 1)
group2 = which(group == 2)
meta2 = cbind(meta1, group)

ttest_ns <- function(exp, group1, group2) {
  tstat = t.test(exp[group1], exp[group2], alternative="less")
  as.numeric(tstat[1:3]) 
}
t01 = t(apply(d2, 1, ttest_ns, group1=group1, group2=group2))
p_none = t01[,3]
p_bonf = p.adjust(p_none, method="bonferroni")
p_bh = p.adjust(p_none, method="BH")
p = cbind(not_adjusted=p_none, Bonferroni_adjusted=p_bonf, BH_FDR_adjusted=p_bh)
sig.genes = p <= 0.05

library(limma)
vennDiagram(vennCounts(sig.genes), include="both", cex = 1, counts.col='red')

library(siggenes)
sam.out = sam(d2, group, rand=123)

print(sam.out, seq(0.1, 5.0, 0.1))

findDelta(sam.out, fdr = 0.05)
plot(sam.out, 1.323008)

findDelta(sam.out, genes = 7000)
plot(sam.out, 4.009168)

sum.sam.out = summary(sam.out, 4.009168)
#print(sum.sam.out, varNames="Probes")
sam.sig.genes = sum.sam.out@row.sig.genes

sig.genes = cbind(sig.genes, SAM=0)
sig.genes[sam.sig.genes, 4] = 1
vennDiagram(vennCounts(sig.genes[,c(2,3,4)]), include="both", names=c("Bonferroni_adjusted", "BH_FDR_adjusted", "SAM[Delta=4.009168]"), cex = 1, counts.col='red')


d2_a = d2[which(sig.genes[,2] == 1), ]
d2_b = d2[p[,2] <= 0.0005, ]

#scale
d3 <- t(scale(t(d2_a), center=TRUE, scale=TRUE))
wss = (nrow(d3) -1) * sum(apply(d3, 2, var))
for (i in 2:15) {
  wss[i] = sum(kmeans(d3, centers=i)$withinss)
}
plot(1:15, wss, type='b', xlab="# of clusters [k]", ylab="Within groups sum of squares")

cl1.1.group = kmeans(d3, 10)$cluster

cordist <- function(x) as.dist(1-cor(t(x), method="pearson"));

cl2.1 = hclust(dist(d3, method="euclidean"), method='ward');
cl2.1.group <- cutree(cl2.1, k=10);
cl2.2 = hclust(cordist(d3), method='ward');
cl2.2.group <- cutree(cl2.2, k=10);

library(pvclust)
cl3 <- pvclust(d3, method.hclust="ward", method.dist="correlation")
plot(cl3) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(cl3, alpha=.95) 

library(mclust)
cl3 <- Mclust(d3)
plot(cl3, d3) # plot results
print(cl3) # display the best model 

adjustedRandIndex(cl2.1.group, cl2.2.group)

library(cluster)
clusplot(d3, cl1.1$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)

library(fpc)
plotcluster(d3, cl1.1.group, xlab='discriminant coordinate 1', ylab='discriminant coordinate 2') 
cluster.stats(dist(d3, method="euclidean"), cl1.1.group, cl1.2.group) 

et = data.frame(as.table(d3))
colnames(et) = c("probe", "tissue", "expression")

hc = cl2.2

mapping1 = data.frame(idx=hc$order, treeIdx=1:nrow(d3))
mapping2 = cbind(mapping1[order(mapping1$idx),], probe=rownames(d3))
et2 = merge(et, mapping2, by.x="probe", by.y="probe")
p <- ggplot(et2, aes(treeIdx, tissue)) + geom_tile(aes(fill=expression)) + 
  scale_fill_gradient(low='white', high='red') +
  scale_y_discrete(name='hybridization', breaks=paste("ME1_H", seq(1,52,3), sep="")) +
  opts(axis.text.x = theme_text(hjust=1, size=10)) +
  opts(axis.text.y = theme_text(hjust=1, size=10));
ggsave(p, filename = file.path(dirO, "project.png"), width=10, height=5);

