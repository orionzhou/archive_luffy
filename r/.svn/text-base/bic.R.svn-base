require(BICseq)
dir_bam = file.path(DIR_Repo, "mt_35/11_pipe_bwa/13_dup_removed")

chrs = paste("chr", 1:8, sep="")
acc1 = "HM001"
accR = "HM101"

f_bam1 = paste(dir_bam, "/", acc1, ".bam", sep="")
f_bamR = paste(dir_bam, "/", accR, ".bam", sep="")

bicseq <- BICseq(sample = f_bam1, reference = f_bamR, seqNames = chrs)
segs <- getBICseg(object = bicseq, bin = 100, lambda = 2, winSize = 200, quant = 0.95, mult = 1)

plot(segs, sampleName = "Demo", save = FALSE, plotBin=TRUE)

bins <- BICseq:::getRatios(bin(segs), what = "bin")
seg.summary <- BICseq:::getSummary(segs, correction=TRUE)



