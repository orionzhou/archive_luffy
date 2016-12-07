
pid = "grn23"
dirw = file.path("/home/springer/zhoux379/scratch/shortread", pid)

fi = file.path(dirw, '00.5.htseq.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

gids = c()
for (i in 1:nrow(ti)) {
	fh = ti$HtseqFile[i]
	th = read.table(fh, header = F, sep = "\t", as.is = T)
	ngene = nrow(th) - 5
	gids0 = th$V1[-c(ngene+1:5)]
	vals = th$V2[-c(ngene+1:5)]
	vals = vals / (sum(vals) / 1000000)
	if(i == 1) {
		gids = gids0
		to = data.frame(gid = gids, vals = vals, stringsAsFactors = F)
	} else {
		stopifnot(identical(gids, gids0))
		to = cbind(to, vals = vals)
	}
	colnames(to)[i+1] = ti$SampleID[i]
}

dim(to)
fo = file.path(dirw, '33.htseq.tsv')
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
