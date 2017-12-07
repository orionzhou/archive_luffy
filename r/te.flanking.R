require(plyr)
require(ggplot2)
require(dplyr)

dirw = '/home/springer/zhoux379/data/misc1/te.flanking'

# cd $misc1/te.flanking
# seq.splitn.py --n 8 --cpu 8 02.fas 04.blat.W22
# seq.splitn.py --n 8 --cpu 8 02.fas 04.blat.PH207
# cat 04.blat.W22/*.psl > 05.blat.W22.psl
# cat 04.blat.PH207/*.psl > 05.blat.PH207.psl
# psl2tsv.pl -i 05.blat.W22.psl -o 06.W22.tsv
# psl2tsv.pl -i 05.blat.PH207.psl -o 06.PH207.tsv
# psl.filter.py --ident 0.9 --cov 0.9 06.W22.tsv 07.W22.tsv
# psl.filter.py --ident 0.9 --cov 0.9 06.PH207.tsv 07.PH207.tsv

fl = file.path(dirw, "02.tsv")
tl = read.table(fl, sep = "\t", as.is = T, header = F)
colnames(tl) = c("id", "size")
ids_all = tl$id

fi = file.path(dirw, "07.W22.tsv")
tir = read.table(fi, sep = "\t", header = T, as.is = T)[,1:21]
ti = tir[tir$alnLen/tir$qSize >= 0.9 & tir$ident >= 0.9,]

ids_mapped = unique(ti$qId)

grp = dplyr::group_by(ti, qId)
ti2 = dplyr::summarise(grp, nbest = sum(score == max(score)), score = max(score))
ti3 = merge(ti, ti2, by = c("qId", "score"))
ids_mapped_u = unique(ti2$qId[ti2$nbest == 1])
ids_mapped_m = unique(ti2$qId[ti2$nbest > 1])

# unique hits
tk = ti3[ti3$qId %in% ids_mapped_u,]
stopifnot(nrow(tk) == length(ids_mapped_u))
ids_mapped_u1 = tk$qId[tk$alnLen == tk$qSize & tk$misMatch == 0 & tk$qNumIns+tk$tNumIns==0]
ids_mapped_u2 = tk$qId[tk$alnLen < tk$qSize | tk$misMatch > 0 | tk$qNumIns+tk$tNumIns > 0]
stopifnot(length(ids_mapped_u1) + length(ids_mapped_u2) == length(ids_mapped_u))


outputs = c(
'',
sprintf("%5d total sequences:", length(ids_all)),
sprintf("%5d do not map to W22 (0.9 identity, 0.9 coverage)", length(ids_all)-length(ids_mapped)),
sprintf("%5d mapped >=1 times:", length(ids_mapped)),
sprintf("   %5d mapped uniquely:", length(ids_mapped_u)),
sprintf("      %5d are identical:", length(ids_mapped_u1)),
sprintf("      %5d have at least 1 mismatch/indel:", length(ids_mapped_u2)),
sprintf("   %5d mapped multiple times:", length(ids_mapped_m)),
''
)
cat(paste(outputs, collapse = "\n"))
