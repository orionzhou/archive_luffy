library(plyr)
library(rtracklayer)
library(Cairo)
library(GenomicRanges)
library(VennDiagram)

dir = '/home/youngn/zhoup/Data/misc3/pan4seq'

# venn-diagram of novel sequences shared
fc = file.path(dir, '41.tbl')
tc = read.table(fc, header = T, sep = "\t", as.is = T)
colnames(tc) = c("cid", "len", 'src', 'org', "orgs", "cnts", "strs")

tc$src[which(tc$src == 'unc')] = 'plant'
ddply(tc, .(orgs), summarise, total_len = sum(len))
x = ddply(tc[tc$len >= 50,], .(orgs), summarise, total_len = sum(len))
ddply(tc[tc$len >= 1000,], .(orgs), summarise, total_len = sum(len))
y = ddply(tc[tc$len >= 50,], .(orgs, src), summarise, total_len = sum(len))


labels = c("HM340", "HM034", "HM056")
n1 = x$total_len[x$orgs == 'HM340.APECCA']
n2 = x$total_len[x$orgs == 'HM034']
n3 = x$total_len[x$orgs == 'HM056']
n12 = x$total_len[x$orgs == 'HM340.APECCA,HM034']
n13 = x$total_len[x$orgs == 'HM340.APECCA,HM056']
n23 = x$total_len[x$orgs == 'HM034,HM056']
n123 = x$total_len[x$orgs == 'HM340.APECCA,HM034,HM056']

venn.plot <- draw.triple.venn(
  area1 = n1+n12+n13+n123, 
  area2 = n2+n12+n23+n123,
  area3 = n3+n13+n23+n123,
  n12 = n12+n123,
  n23 = n23+n123,
  n13 = n13+n123,
  n123 = n123,
  category = labels,
  fill = c("blue", "red", "limegreen"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("blue", "red", "limegreen"),
  alpha = 0.5,
  euler.d = T,
  scaled = T,)
CairoPNG(filename = file.path(dir, "venn.png"), width = 600, height = 600)
grid.draw(venn.plot)
dev.off()

# plot novel segments length distribution
tmp1 = table(tc$len)
tp.1 = data.frame(len=as.numeric(names(tmp1)), cnt=c(tmp1))
tp.2 = cbind(tp.1, sum=tp.1$len * tp.1$cnt)
tp.3 = tp.2[order(tp.2$len, decreasing=T),]
tp = cbind(tp.3, cumsum=cumsum(tp.3$sum))
plot(tp$len, tp$cumsum, type='l', xlim=c(0, 8000), main='novel', xlab='segment length', ylab='cumsum of segments')
x=c(50,100,1000)
y=c(tp$cumsum[tp$len==50], tp$cumsum[tp$len==100], tp$cumsum[tp$len==1000])
segments(0, y, x, y, col='blue')
segments(x, y, x, 0, col='blue')


##### blast NR/NT database
fb = file.path(dir, "31_blastnr/15.gal")
tbl = read.table(fb, header = T, sep = "\t", as.is = T)[,c(1:12,18,21:25)]
tb = ddply(tbl, .(qId), summarise, 
  ali=sum(ali), n_seg=length(tId), score = sum(score),
  n_tax=length(unique(cat)),
  cat=cat[which(score==max(score))[1]], species=species[1])

t41 = cbind(t21c, category='', stringsAsFactors=FALSE)
t41$category[t21c$qLen>0] <- 'Medicago'
for (i in 1:nrow(t37)) {
  j = which(t41$qId == t37$qId[i])
  if(t41$score[j] < t37$score[i]) {
    t41[j, 2:7] = t37[i, 2:7]
    t41$category[j] = t37$category[i]
    t41$pct_cov[j] = t41$qLen[j] / t41$len_scaf[j]
  }
}
table(t41$category)

scaffold_stat <- function(ids, df) {
  dfs = df[df$qId %in% ids,]
  c('n'=length(ids), 'len_scaf'=sum(dfs$len_scaf), 'qLen_aln'=sum(dfs$qLen), 'hLen_aln'=sum(dfs$hLen))
}
id_all = t41$qId
id_mt = t41$qId[t41$category == "Medicago"]
id_fa = t41$qId[t41$category == "Fabaceae"]
id_un = t41$qId[t41$category == ""]
id_misc = id_all[! id_all %in% c(id_mt, id_fa, id_un)]
sapply(list(all=id_all, medicago=id_mt, fabaceae=id_fa, unknown=id_un, misc=id_misc), scaffold_stat, t41)

f41 = file.path(dir, "41_scaffold_status.tbl")
write.table(t41, file=f41, col.names=T, row.names=F, sep="\t", quote=F)
t41b = t41[t41$qId %in% c(id_fa, id_un),]
write.table(t41b, file=file.path(dir, "42_unknown.tbl"), col.names=T, row.names=F, sep="\t", quote=F)



# construct pseudo-chrs
# get scaffold order
tm = read.table(file.path(dir, "25_blat_final/35.gal"), header=T, sep="\t", as.is=T)[,1:17]
get_row_max_score <- function(df) { df[which.max(df[,'score']),] }
tm.s = ddply(tm, .(tId), get_row_max_score)
to.1 = data.frame(id=tm.s$tId, chr=tm.s$qId, pos=(tm.s$qBeg+tm.s$qEnd)/2, srd=tm.s$qSrd, idx=0, stringsAsFactors=F)
to = merge(to.1, ts.q, by='id', all=T)[,-6]
to = to[order(to$chr, to$pos, to$id),]
to$idx = 1:nrow(to)

tt = data.frame(id=to$id, cat_blat='mt', idx=to$idx, stringsAsFactors=F)
tt$cat_blat[is.na(to$chr)] = rep('unc', sum(is.na(to$chr)))

# blast NT
tb = read.table(file.path(dir, "51_pan3/15.tbl"), header=T, sep="\t")[,-c(18,19)]
tb.1 = ddply(tb, .(qId, cat), summarise, score=sum(score))
tb.2 = ddply(tb.1, .(qId), get_row_max_score)
tb = data.frame(id=tb.2$qId, cat_nt=as.character(tb.2$cat), stringsAsFactors=F)

# classify novelseq
ti = read.table(file.path(dir, "51_pan3/01.tbl"), header=T, sep="\t")[,-c(18,19)]
colnames(ti)[1] = 'chr'
makelocid <- function(x) sprintf("%s-%d-%d", x['chr'], as.numeric(x['beg']), as.numeric(x['end']))
ti = cbind(id=apply(ti, 1, makelocid), ti)

ti.1 = merge(ti, tb, by="id", all=T)
ti.2 = merge(ti.1, tt, by.x='chr', by.y='id', all.x=T)
ti.2$cat_nt[which(is.na(ti.2$cat_nt))] <- 'unc'
table(ti.2[,c('cat_nt', 'cat_blat')])
ti.3 = cbind(ti.2, cat='plant', stringsAsFactors=F)
ti.3$cat[which(ti.2$cat_blat=='unc' & ti.2$cat_nt=='foreign')] = 'unc'

tf1 = ti.3[ti.3$cat=='plant',]
tf1 = tf1[order(tf1$idx), c(1,3,4,5)]
cat(org.q, "specific:", nrow(tf1), "segments,", sum(tf1$len), "bp\n")
write.table(tf1, file.path(dir, "51_pan3/21_plant.tbl"), col.names=T, row.names=F, sep='\t', quote=F)
tf2 = ti.3[ti.3$cat=='unc',]
tf2 = tf2[order(tf2$idx), c(1,3,4,5)]
cat(org.q, "unc:", nrow(tf2), "segments,", sum(tf2$len), "bp\n")
write.table(tf2, file.path(dir, "51_pan3/22_unc.tbl"), col.names=T, row.names=F, sep='\t', quote=F)
