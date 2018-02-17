require(tidyverse)

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

### batch 1 sample mix-up correction
fi = file.path(dirw, "10.read.tsv")
ti = read.table(fi, sep = "\t", header = T, stringsAsFactors = F)

sms1 = c("BR003", "BR006", "BR032")
sms2 = c("BR004", "BR007", "BR029")

for (i in 1:length(sms1)) {
    sm1 = sms1[i]; sm2 = sms2[i]
    idx1 = which(ti$SampleID == sm1)
    idx2 = which(ti$SampleID == sm2)
    gt1 = ti$Genotype[idx1]; gt2 = ti$Genotype[idx2]
    ti$Genotype[idx1] = gt2; ti$Genotype[idx2] = gt1
}

to = ti
to$Tissue = factor(to$Tissue, levels = unique(to$Tissue))
to$Genotype = factor(to$Genotype, levels = unique(to$Genotype))
to = to[order(to$Tissue, to$Genotype, to$SampleID),]

tos = unique(to[,c("Tissue",'Genotype')])
for (i in 1:nrow(tos)) {
    tiss = tos$Tissue[i]; gt = tos$Genotype[i]
    idxs = which(to$Tissue == tiss & to$Genotype == gt)
    to$Treatment[idxs] = 1:length(idxs)
}
fo = file.path(dirw, "11.read.correct.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)


### batch 2 sample identification
fi = file.path(dirw, "11.read.correct.tsv")
ti = read.table(fi, sep = "\t", header = T, stringsAsFactors = F)[,1:5]

fd = file.path(dirw, "exp.design.tsv")
td = read.table(fd, sep = "\t", header = F, stringsAsFactors = F)
colnames(td) = c("SampleID", "Tissue", "Genotype", "Treatment", "Pool", "Index")
td = cbind(td, cood = sprintf("%s_%s", tolower(td$Pool), substr(td$Index, nchar(td$Index)-1, nchar(td$Index))))

f_gt = file.path(dirw, "../46.ase", "13.gt.tsv")
t_gt = read_tsv(f_gt)

ti2 = merge(ti, td[,c("SampleID","Pool","Index")], by = 'SampleID')
ti3 = merge(ti2, t_gt, by = 'SampleID')
ti4 = ti3 %>%
    filter(SampleID > 'BR165') %>%
    arrange(Tissue, Genotype, Treatment)


tdic = c(
    "BR197" = 'seedlingmeristem_11DAS',
    "BR198" = 'seedlingmeristem_11DAS',
    'BR215' = "seedlingroot_11DAS",
    'BR216' = 'seedlingroot_11DAS',
    'BR187' = 'radicle_root',
    'BR227' = 'radicle_root',
    'BR212' = 'radicle_root',
    'BR214' = 'radicle_root',
    'BR179' = 'seedlingleaf_11DAS',
    'BR180' = 'seedlingleaf_11DAS'
)
df_tdic = data.frame(SampleID = names(tdic), tiss = as.character(tdic), stringsAsFactors = F)
ti5 = merge(ti4, df_tdic, by = 'SampleID', all.x = T)
ti5$tiss[is.na(ti5$tiss)] = ti5$Tissue[is.na(ti5$tiss)]
ti5 = ti5 %>% 
    mutate(tiss_ok = ifelse(Tissue == tiss, '1', '3')) %>%
    arrange(tiss, gt)

fo = file.path(dirw, "15.batch2.tsv")
write.table(ti5, fo, sep = "\t", row.names = F, col.names = T, quote = F)


### batch 2 sample label correction
fi = file.path(dirw, "11.read.correct.tsv")
ti = read.table(fi, sep = "\t", header = T, stringsAsFactors = F)

f_gt = file.path(dirw, "../46.ase", "13.gt.tsv")
t_gt = read_tsv(f_gt)

tdic = c(
    "BR197" = 'seedlingmeristem_11DAS',
    "BR198" = 'seedlingmeristem_11DAS',
    'BR215' = "seedlingroot_11DAS",
    'BR216' = 'seedlingroot_11DAS',
    'BR187' = 'radicle_root',
    'BR227' = 'radicle_root',
    'BR212' = 'radicle_root',
    'BR214' = 'radicle_root',
    'BR179' = 'seedlingleaf_11DAS',
    'BR180' = 'seedlingleaf_11DAS'
)
df_tdic = data.frame(SampleID = names(tdic), tiss = as.character(tdic), stringsAsFactors = F)

ti2 = merge(ti, t_gt, by = 'SampleID')
ti3 = merge(ti2, df_tdic, by = 'SampleID', all.x = T)
ti3$tiss[is.na(ti3$tiss)] = ti3$Tissue[is.na(ti3$tiss)]

to = transmute(ti3, SampleID=SampleID, Species=Species, Tissue=tiss, Genotype=gt, 
    Treatment=Treatment, ReadFile1=ReadFile1, ReadFile2=ReadFile2)
to$Genotype = factor(to$Genotype, levels = unique(to$Genotype))
to$Tissue = factor(to$Tissue, levels = unique(ti3$Tissue))
to = arrange(to, Tissue, Genotype, SampleID)

tos = unique(to[,c("Tissue",'Genotype')])
for (i in 1:nrow(tos)) {
    tiss = tos$Tissue[i]; gt = tos$Genotype[i]
    idxs = which(to$Tissue == tiss & to$Genotype == gt)
    to$Treatment[idxs] = 1:length(idxs)
}


fo = file.path(dirw, "12.read.batch2correct.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
