require(plyr)
require(dplyr)
require(tidyr)
require(ggplot2)

dirw = '/home/springer/zhoux379/scratch/briggs'
dirw = '/home/springer/zhoux379/data/misc2/briggs/01.exp.design'

fi = file.path(dirw, "readfiles.txt")
ti = read.table(fi, stringsAsFactors = F)
x = strsplit(ti$V1, split = "_")
pool = sapply(x, "[", 2)
index = sapply(x, "[", 4)
pair = sapply(x, "[", 6)
tf = data.frame(fn = ti$V1, pair = pair, cood = sprintf("%s_%02d", pool, as.integer(substr(index, 6, nchar(index)))), stringsAsFactors = F)
tf2 = spread(tf, pair, fn)

fd = file.path(dirw, "exp.design.tsv")
td = read.table(fd, sep = "\t", header = F, stringsAsFactors = F)
colnames(td) = c("SampleID", "Tissue", "Genotype", "Treatment", "Pool", "Index")
td = cbind(td, cood = sprintf("%s_%s", tolower(td$Pool), substr(td$Index, nchar(td$Index)-1, nchar(td$Index))))

tm = merge(td, tf2, by = 'cood', all.x = T)
tm = tm[order(tm$SampleID),]
to = data.frame(SampleID = tm$SampleID, Species = "Zmays", Tissue = tm$Tissue, Genotype = tm$Genotype, Treatment = tm$Treatment, ReadFile1 = tm$R1, ReadFile2 = tm$R2, stringsAsFactors = F)

#to = to[!is.na(to$ReadFile1),]
for (i in 1:nrow(to)) {
	f1a = file.path(dirw, "01.reads.1", to$ReadFile1[i])
	f1b = file.path(dirw, "01.reads.2", to$ReadFile1[i])
	f1c = file.path(dirw, "01.reads.3", to$ReadFile1[i])
	f1d = file.path(dirw, "01.reads.4", to$ReadFile1[i])
	f2a = file.path(dirw, "01.reads.1", to$ReadFile2[i])
	f2b = file.path(dirw, "01.reads.2", to$ReadFile2[i])
	f2c = file.path(dirw, "01.reads.3", to$ReadFile2[i])
	f2d = file.path(dirw, "01.reads.4", to$ReadFile2[i])
	if(file.exists(f1a) & file.exists(f2a)) {
		to$ReadFile1[i] = normalizePath(f1a)
		to$ReadFile2[i] = normalizePath(f2a)
	} else if (file.exists(f1b) & file.exists(f2b)) {
		to$ReadFile1[i] = normalizePath(f1b)
		to$ReadFile2[i] = normalizePath(f2b)
	} else if (file.exists(f1c) & file.exists(f2c)) {
		to$ReadFile1[i] = normalizePath(f1c)
		to$ReadFile2[i] = normalizePath(f2c)
	} else if (file.exists(f1d) & file.exists(f2d)) {
		to$ReadFile1[i] = normalizePath(f1d)
		to$ReadFile2[i] = normalizePath(f2d)
	} else {
	    cat(to$SampleID[i], to$ReadFile1[i], "\n")
	}
}
fo = file.path(dirw, "10.read.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### sample mix-up correction
fi = file.path(dirw, "10.read.tsv")
ti = read.table(fi, sep = "\t", header = T, stringsAsFactors = F)

sms1 = c("BR003", "BR006", "BR032")
sms2 = c("BR004", "BR007", "BR029")

for (i in 1:length(sms1)) {
    sm1 = sms1[i]; sm2 = sms2[i]
    idx1 = which(ti$SampleID == sm1)
    idx2 = which(ti$SampleID == sm2)
    ti$SampleID[idx1] = sm2
    ti$SampleID[idx2] = sm1
    f1 = ti$ReadFile1[idx1]; f2 = ti$ReadFile1[idx2]
    ti$ReadFile1[idx1] = f2; ti$ReadFile1[idx2] = f1
    f1 = ti$ReadFile2[idx1]; f2 = ti$ReadFile2[idx2]
    ti$ReadFile2[idx1] = f2; ti$ReadFile2[idx2] = f1
}
fo = file.path(dirw, "11.read.correct.tsv")
write.table(ti, fo, sep = "\t", row.names = F, col.names = T, quote = F)
