require(seqinr)

tname = "HM101"
qnames = c(
  "HM058", "HM125", "HM056.AC", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340.AC"
)
diro = file.path(Sys.getenv("misc3"), "comp.ortho")

fi = file.path(diro, "33.ortho.cat.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti[is.na(ti)] = ''

orgs = c(tname, qnames)
n_org = apply(ti, 1, function(z) sum(z[1:length(orgs)] != ''))
ti = cbind(ti, n_org = n_org)


ids = c('Medtr4g020620.1', 'Medtr8g060730.1', 'Medtr5g088410.1', 'Medtr2g020630.1', 'Medtr2g044570.1')

to = ti[ti[,tname] %in% ids,]
rownames(to) = to[,tname]

fass = list()
for (org in c(tname, qnames)) {
  f_fas = file.path(Sys.getenv("genome"), org, "51.fas")
  fas <- read.fasta(f_fas, seqtype = "AA", as.string = T, set.attributes = F)
  fass[[org]] = fas
}


for (id in ids) {
  seqs = list()
  for (org in c(tname, qnames)) {
    nid = to[id, org]
    if(id != '') {
      seqs[paste(org, nid, sep = "_")] = fass[[org]][[nid]]
    }
  }
  fos = sprintf("%s/cases/%s.fas", diro, id)
  write.fasta(sequences = seqs, names = nids, nbchar = 80, file.out = fos)
}