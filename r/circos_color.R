fi = file.path(DIR_Misc1, "circos/01_conf/circos_color.tbl")
d1 = read.table(fi, header=F, stringsAsFactors=F, sep="\t")
colnames(d1)=c('col', 'num', 'type', 'cnt', 'r', 'g', 'b')
d2 = d1[!is.na(as.numeric(d1$cnt)),]
d2$num = as.numeric(d2$num)
d2$cnt = as.numeric(d2$cnt)

getLabel <- function(x) { paste(x[c(2,4,3)], sep="-", collapse="-") }
p1 = unique(d2[,c('col','num','type')])
p2 = aggregate(p1$num, by=list(factor(p1$col), factor(p1$type)), FUN=max)
colnames(p2) = c('col', 'type', 'num')
p3 = cbind( pal=1:nrow(p2), p2[order(p2$col,p2$num),] )
rownames(p3) = 1:nrow(p3)
p4 = cbind( p3, label=apply(p3, 1, getLabel) )

d4 = merge(p4, d2, by=c('col', 'num', 'type'))
d4 = d4[order(d4$pal, d4$cnt),]
rownames(d4) =1:nrow(d4)
getRGB <- function(x) {rgb(x[6], x[7], x[8], maxColorValue=255)}
d5 = cbind(d4, rgb_cnt=1:nrow(d4), rgb = apply(d4, 1, getRGB))

png(file=file.path(DIR_Misc1, "circos/01_conf/circos_color.png"), width=900, height=900)
plot(c(0, max(d5$cnt)), c(0, max(d5$pal)), type='n', xlab='seq', ylab='', yaxt='n', main='')
i = 1:nrow(d5)
rect(d5$cnt[i]-1, d5$pal[i]-1, d5$cnt[i], d5$pal[i], col=rgb(d5$r[i],d5$g[i],d5$b[i],maxColorValue=255), border=NA)
axis(2, at=p4$pal-0.5, labels=p4$label, col.axis="blue", las=2, pos=0, cex.axis=0.9)
dev.off()


