library(Rgraphviz)

org = "Athaliana"
dir = file.path(DIR_Misc2, "crp.evo", org)
f01 = file.path(dir, "02_seq.tbl")
ts = read.table(f01, header=T, sep="\t", as.is=T)

f05 = file.path(dir, "05_pairs.tbl")
tp = read.table(f05, header=T, sep="\t", as.is=T)

colnames(tp) = c("from", "to", "weight")
g = graphBAM(tp)
cc = connComp(g)

clus = rep(sprintf("clu%03d",1:length(cc)), unlist(lapply(cc,length)))
cd = data.frame(id=unlist(cc), clu=clus)
dc = merge(ts[,1:7], cd, by="id")

