f01 = file.path("/project/youngn/zhoup/Data/misc3/gm_snp1/goldengate1536.csv")
d1 = read.table(f01, sep=",", header=TRUE)

f02 = file.path("/project/youngn/zhoup/Data/misc3/gm_snp1/04_snp_loc.txt")
d2 = read.table(f02, sep="\t", header=FALSE)
colnames(d2) = c("id", "allele", "chr", "pos", "pct_idty", "alias", "note")

d3 = merge(d1, d2, by.x="Locus", by.y="id", all.x=TRUE, all.y=FALSE)
write.table(d3, file.path("/project/youngn/zhoup/Data/misc3/gm_snp1/goldengate1536_location.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


