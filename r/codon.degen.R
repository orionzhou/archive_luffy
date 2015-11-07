require(seqinr)

chs = c("A", "C", 'G', "T")
for (i in chs) {
  for (j in chs) {
    for (k in chs) {
      cod = paste(i, j, k, sep = '')
      aln = as.alignment(nb = 2, nam = c('s1', 's2'), seq = c(cod, cod))
      res = kaks(aln, verbose = T)
      l0 = res$l0[1]; l2 = res$l2[1]; l4 = res$l4[1]
      cat(cod, l0, l2, l4, "\n", sep = " ")
    }
  }
}