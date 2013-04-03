library(data.table)

#acc = "hm056"
acc = "hm340"
dir = file.path(DIR_Misc3, acc)

f_len = file.path(dir, "11_seqlen.tbl")
t_len = read.table(f_len, header=TRUE, sep="\t", as.is=T)
colnames(t_len)=c("qId", "len_scaf")

f_gap = file.path(dir, "12_gaploc.tbl")
t_gap = read.table(f_gap, header=TRUE, sep="\t", as.is=T)
colnames(t_gap)=c("qId", "qBeg", "qEnd", "len")

# calculate assembly statistics
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/contigStats.R") 
library(Biostrings)
assembly <- readDNAStringSet(sprintf("%s/01_assembly.fa", dir, acc), "fasta")
N <- list(acc=width(assembly))
reflength <- sapply(N, sum)
stats <- contigStats(N=N, reflength=reflength, style="data")
stats[["Contig_Stats"]]

dir = file.path(dir, '21_blastn')

f05 = file.path(dir, "05_tiled.tbl")
t05 = read.table(f05, header=TRUE, sep="\t", as.is=T)
t05  = t05[order(t05$qId, t05$qBeg),]

# get unmapped scaffolds for NR blast
t21 = ddply(t05, .(qId), summarise, qLen=sum(qLen), hLen=sum(hLen), n_seg=length(hId), pct=mean(pct), e=min(e), score = sum(score))
t21b = merge(t21, t_len, by="qId", all=T)
t21b[is.na(t21b)] = 0
t21c = cbind(t21b, pct_cov=t21b$qLen/t21b$len_scaf)
t21d = t21c[t21c$qLen<1000,]

f21 = file.path(dir, "21_scaffold_status.tbl")
write.table(t21c, file=f21, col.names=T, row.names=F, sep="\t", quote=F)
write.table(t21d, file=file.path(dir, "22_unmapped.tbl"), col.names=T, row.names=F, sep="\t", quote=F)

# run assembly.pl to blast NR/NT database

f36 = file.path(dir, "36_annotated.tbl")
t36 = read.table(f36, header=TRUE, sep="\t", as.is=T)
t37 = ddply(t36, .(qId), summarise, qLen=sum(qLen), hLen=sum(hLen), n_seg=length(hId), pct=mean(pct), e=min(e), score = sum(score), n_tax=length(unique(category)), category=category[which(score==max(score))[1]], species=species[1])

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


