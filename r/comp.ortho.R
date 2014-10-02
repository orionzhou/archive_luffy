require(rtracklayer)
require(plyr)

args <- commandArgs(trailingOnly = TRUE)
print(args)

qname <- ifelse(is.na(args[1]), "HM004", as.character(args[1]))
tname = "HM101"

dir = sprintf("%s/%s_%s/23_blat", Sys.getenv("misc3"), 
  toupper(qname), toupper(tname))
fn = file.path(dir, "ortho.tbl")
tn = read.table(fn, sep = "\t", header = T, as.is = T)
tns = tn[tn$qid != '' & tn$tid != '',]

aa_pw_dist <- function(rw) {
s1 = rw['tseq']
s2 = rw['qseq']
pw = pairwiseAlignment(AAString(s1), AAString(s2), substitutionMatrix = "BLOSUM62", gapOpening = -3, gapExtension = -1)
len = nchar(pw)
indel = nindel(pw)
mis = nrow(mismatchTable(pw)) + indel@insertion[2] + indel@deletion[2]
mis / len
}

cl = makeCluster(detectCores())

cluster_fun <- function() {
    require(Biostrings)
    require(plyr)
}
clusterCall(cl, cluster_fun)


ptm <- proc.time()
y = parApply(cl, tns, 1, aa_pw_dist)
proc.time() - ptm


stopCluster(cl)

to = cbind(tns[,1:5], ident = y)
fo = file.path(dir, "ortho.2.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
