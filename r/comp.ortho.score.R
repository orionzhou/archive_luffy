require(rtracklayer)
require(plyr)
require(seqinr)
source("Align.R")
source("comp.fun.R")

#args <- commandArgs(trailingOnly = TRUE)
#print(args)
#qname <- ifelse(is.na(args[1]), "HM004", as.character(args[1]))

##### cluster initialization
cl = makeCluster(detectCores())

cluster_fun <- function() {
    require(Biostrings)
    require(plyr)
}
clusterCall(cl, cluster_fun)

##### determine ortholog identity using pairwise alignment
tname = "HM101"
qnames = get_orgs()

for (qname in qnames) {
#qname = "HM004"
f_tfas = file.path(Sys.getenv("genome"), tname, "51.fas")
tfas <- read.fasta(f_tfas, seqtype = "AA", as.string = T, set.attributes = F)
f_qfas = file.path(Sys.getenv("genome"), qname, "51.fas")
qfas <- read.fasta(f_qfas, seqtype = "AA", as.string = T, set.attributes = F)

diro = sprintf("%s/%s_%s/51_ortho", Sys.getenv("misc3"), qname, tname)

fi = file.path(diro, "01.ortho.tbl")
ti = read.table(fi, sep = "\t", header = T, as.is = T)[,1:4]

tm = cbind(ti, tseq = as.character(tfas[ti$tid]), 
  qseq = as.character(qfas[ti$qid]), stringsAsFactors = F)

ptm <- proc.time()
y = parApply(cl, tm, 1, aa_pw_dist)
proc.time() - ptm

to = cbind(ti, t(y))
to$qlen = as.integer(to$qlen / 3)
to$tlen = as.integer(to$tlen / 3)
fo = file.path(diro, "05.score.tbl")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
}

##### stop cluster
stopCluster(cl)