require(plyr)
require(dplyr)
require(ggplot2)

dirw = '/home/springer/zhoux379/scratch/briggs'

fi = file.path(dirw, "filenames.txt")
ti = read.table(fi, stringsAsFactors = F)
x = strsplit(ti$V1, split = "_")
pool = sapply(x, "[", 2)
index = sapply(x, "[", 4)
pair = sapply(x, "[", 6)
tf = data.frame(fn = ti$V1, pair = pair, cood = sprintf("%s_%02d", pool, as.integer(substr(index, 6, nchar(index)))), stringsAsFactors = F)
tf2 = reshape(tf, direction = "wide", idvar=c("cood"), timevar="pair")

fd = file.path(dirw, "exp.design.tsv")
td = read.table(fd, sep = "\t", header = T, stringsAsFactors = F)
td = cbind(td, cood = sprintf("%s_%s", tolower(td$Pool), substr(td$Index, nchar(td$Index)-1, nchar(td$Index))))

tm = merge(td, tf2, by = 'cood', all.x = T)
tm = tm[order(tm$ID),]
to = data.frame(SampleID = tm$Label, Species = "Zmays", Tissue = tm$Tissue, Genotype = tm$Genotype, Treatment = tm$Replicate, ReadFile1 = tm$fn.R1, ReadFile2 = tm$fn.R2, stringsAsFactors = F)

#to = to[!is.na(to$ReadFile1),]
for (i in 1:nrow(to)) {
	f11 = file.path(dirw, "05.reads.1", to$ReadFile1[i])
	f12 = file.path(dirw, "05.reads.1", to$ReadFile2[i])
	f21 = file.path(dirw, "05.reads.2", to$ReadFile1[i])
	f22 = file.path(dirw, "05.reads.2", to$ReadFile2[i])
	if(file.exists(f11) & file.exists(f12)) {
		to$ReadFile1[i] = file.path("05.reads.1", to$ReadFile1[i])
		to$ReadFile2[i] = file.path("05.reads.1", to$ReadFile2[i])
	} else {
		stopifnot(file.exists(f21) & file.exists(f22))
		to$ReadFile1[i] = file.path("05.reads.2", to$ReadFile1[i])
		to$ReadFile2[i] = file.path("05.reads.2", to$ReadFile2[i])
	}
}
fo = file.path(dirw, "00.1.read.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

