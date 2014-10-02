require(ape)
require(geiger)
require(igraph)


#### plot NBS-LRR tree (with label)
dir = file.path(Sys.getenv('misc2'), 'nbs/mt_40')
fi = file.path(dir, "11.phb")
fo = file.path(dir, "11.png")

fm = file.path(dir, "../mt_10/13.flt.gal")
tm = read.table(fm, header = T, sep = "\t", as.is = T, quote = "")[, c(2,7)]

fc = file.path(dir, "../mt_10/clade.tbl")
tc = read.table(fc, header = T, sep = "\t", as.is = T, quote = "")[, c(1,4,6:7)]

tmc = merge(tm, tc, by.x = 'qId', by.y = 'tid')

tree = read.tree(fi)

pick_clade <- function(id, df) {
  idx = which(df$tId == id)
  ifelse(length(idx) == 0, '', df$clade[idx[1]])
}
labels = tree$tip.label
clades = sapply(labels, pick_clade, tmc)
tree$tip.label = paste(labels, clades, sep = " ")

  scores = as.numeric(tree$node.label)
  if(mean(scores, na.rm = TRUE) > 1) { scores = scores / 1000 }
  node.labels.bg = rep('white', tree$Nnode)
  node.labels.bg[scores >= 0.95] = 'black'
  node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

png(filename = fo, width = 800, height = 7000, units = 'px')
plot(tree, show.node.label = F, show.tip.label = T, font = 1, 
  label.offset = 0.01, no.margin = T, cex = 0.78)
nodelabels(pch = 22, bg = node.labels.bg)
add.scale.bar(lcol = 'black')
dev.off()

tt = data.frame(id = labels, clade = clades)
write.table(tt, file.path(dir, '15.clade.tbl'), row.names = F, col.names = T,
  sep = "\t", quote = F)

#### plot NBS-LRR tree (simple)
dir = file.path(Sys.getenv('misc2'), 'nbs/mt_40')
fi = file.path(dir, "26.phb")
fo = file.path(dir, "26.png")
tree = read.tree(fi)

png(filename = fo, width = 800, height = 7000, units = 'px')
plot(tree, show.tip.label = T, font = 1, 
  label.offset = 0.01, no.margin = T, cex = 0.78)
add.scale.bar(lcol = 'black')
dev.off()

##### plot crp subgroup tree
source('clustertree.R')

dir = file.path(Sys.getenv("misc2"), 'genefam')
fc = file.path(dir, "11.crp.tbl")
tc = read.table(fc, sep = "\t", header = T, as.is = T)[,1:12]
fo = sprintf("%s/27.stat.tbl", dir)

dat = data.frame()
for (fam in unique(tc$family)) {
  fam = tolower(fam)
  fi = sprintf("%s/24_phb/%s.phb", dir, fam)
  ff = sprintf("%s/26_fig/%s.png", dir, fam)
  if(!file.exists(fi)) next
  
  tree = read.tree(fi)
  cls = prosperi.cluster(tree, 0.05)$membership

  ids = tree$tip.label
  orgs = sapply(strsplit(ids, "_"), "[", 1)

  ucl = sort(unique(cls))
  uorg = sort(unique(orgs))

  pal = rainbow(length(ucl))
  names(pal) = as.character(ucl)

  ht = length(tree$tip.label) * 12
  png(filename = ff, width = 800, height = ht, units = 'px')
  plot(tree, font = 1, show.node.label = F, show.tip.label = F,
    label.offset = 0.01, no.margin = T, cex = 0.81)
  newtip = paste(cls, tree$tip.label, sep = ' ')
  tiplabels(newtip, adj = -0.1, bg = pal[as.character(cls)], cex = 0.8)

  add.scale.bar(lcol = 'black')
  dev.off()

  mat = matrix(0, nrow = length(ucl), ncol = length(uorg), dimnames = 
    list(ucl, uorg))
  for (cl in ucl) {
    tab = table(orgs[cls == cl])
    mat[as.character(cl), names(tab)] = tab
  
    n0 = as.numeric(tab['hm101'])
    n0 = ifelse(is.na(n0), 0, n0)
    for (org in names(tab)) {
      if(org == 'hm101') next
      n = as.numeric(tab[org])
      if((n0 == 0 & n - n0 > 1) | (n0 > 0 & n - n0 > 0))  {
        x = data.frame(fam = fam, cl = cl, hm101 = n0, org = org, n = n,
          ids = paste(ids[cls == cl & orgs == org], collapse = ' '))
        dat = rbind(dat, x)
      }
    }
  }
}
write.table(dat, file = fo, sep = "\t", row.names = F, col.names = T, quote = F)
