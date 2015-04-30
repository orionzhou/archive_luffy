require(GenomicFeatures)

dirw = file.path(Sys.getenv("genome"), "HM101")
fgff = file.path(dirw, "51.gff")

fsize = file.path(dirw, "15.sizes")
fgap = file.path(dirw, "16.gap.bed")
stopifnot(file.exists(fsize) & file.exists(fgap))
  
tlen = read.table(fsize, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))

tgap = read.table(fgap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))

grnp = GenomicRanges::setdiff(grt, grp)

chrominfo = cbind(tlen, V3 = F)
colnames(chrominfo) = c("chrom", "length", "is_circular")

txdb <- makeTranscriptDbFromGFF(file = fgff, format = "gff3",
  exonRankAttributeName = NA, chrominfo = chrominfo, 
  species="Medicago truncatula")

ftxdb = file.path(dirw, "51.sqlite")
if(interactive()) saveDb(txdb, file = ftxdb)