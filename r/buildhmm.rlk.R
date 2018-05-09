require(seqinr)
require(grid)
require(plyr)
require(ggplot2)
require(ape)
source("Align.R")

dirg = file.path(Sys.getenv("genome"), "Athaliana")
dirw = file.path(Sys.getenv("misc2"), 'rlk', 'At')
dir.create(dirw)
setwd(dirw)

fg = file.path(dirg, '51.gtb')
tg = read.table(fg, header = T, as.is = T, sep = "\t")
proseq <- read.fasta(file.path(dirg, "51.fas"), seqtype = "AA", as.string = T, set.attributes = F)

fi = file.path(dirw, "rlk.tbl")
ti = read.table(fi, header = T, as.is = T, sep = "\t")
nrow(ti)

tis = ti[toupper(ti$GeneID) %in% tg$par,]
nrow(tis)
tis = cbind(tis, fam = tis$Subfamily)
tis$fam = gsub(" ", "_", tis$fam)

table(tis$fam)
fams = unique(tis$fam)
fams = fams[! fams %in% c("C-Lectin", "CrRLK1L-2", "N._A.", "N.A.", "RKF3L", "SD-3", "URK_I", "LysM", "Extensin", "Thaumatin")]

tf = tis[tis$fam %in% fams,]
table(tf$fam)

to = cbind(tf, par = toupper(tf$GeneID))
to = merge(to, tg[,c('id','par')], by = 'par')
labs = to$fam
names(labs) = to$id

dir.create('21_fas')
dir.create('22_aln')
dir.create("23_trim")
dir.create("24_hmm")
for (fam in fams) {
  ids = to$id[to$fam == fam]
  pseq = proseq[ids]
  fop = sprintf("21_fas/%s.fas", fam)
  write.fasta(sequences = pseq, names = names(pseq), nbchar = 60, file.out = fop)

  system(sprintf("clustalo -i 21_fas/%s.fas -o 22_aln/%s.aln --outfmt=clu --output-order=tree-order --force --full --full-iter", fam, fam))
  system(sprintf("trimal -automated1 -in 22_aln/%s.aln -out 23_trim/%s.aln", fam, fam))
  system(sprintf("hmmbuild 24_hmm/%s.hmm 23_trim/%s.aln", fam, fam))
#  system(sprintf("clustalw2 -infile=rlk.aln -tree -outputtree=phylip -clustering=NJ -kimura"))
}
system("cat 24_hmm/*.hmm > 30.all.hmm")


### plot a big tree
f_tree = file.path(dirw, "rlk.ph")
tree = read.tree(f_tree)
  ff = sprintf("%s/rlk.png", dirw)
  
  ht = length(tree$tip.label) * 12
  png(filename = ff, width = 1000, height = ht, units = 'px')
  plot(tree, font = 1, show.node.label = F, show.tip.label = F,
    label.offset = 0.01, no.margin = T, cex = 0.81)
  newtip = paste(labs[tree$tip.label], sep = ' ')
  tiplabels(newtip, adj = -0.1, cex = 0.8)

  add.scale.bar(lcol = 'black')
  dev.off()
  