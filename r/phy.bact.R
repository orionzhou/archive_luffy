require(ape)

dirw = file.path(Sys.getenv('misc2'), 'bact')
ft = file.path(dirw, "15.nwk")
tree = read.tree(ft)
tree = root(tree, which(tree$tip.label=='AF075271.2.1482')) #Alterococcus

fl = file.path(dirw, "11.fas.tsv")
tl = read.table(fl, header = F, sep = "\t", as.is = T, quote = "")
colnames(tl) = c("id", "taxa")
labs = sapply(strsplit(tl$taxa, split=";"), "[", 6)
labs = sapply(strsplit(labs, split=" "), "[", 1)
#table(labs)
tl = cbind.data.frame(tl, genus = as.character(labs), stringsAsFactors = F)

  ids = tree$tip.label
  idxs = match(ids, tl$id)
  tl = tl[idxs,]
  
### replace tree label
tree$tip.label = tl$genus
write.tree(tree, file = file.path(dirw, "16.nwk"))

### plot tree
  genus_root = 'Alterococcus'
  genus1 = c('Bifidobacterium', 'Microbispora', 'Actinomadura', 'Kocuria', 'Clavibacter', 'Planomonospora')
  genus2 = c('Actinoplanes')
  genus3 = c('Streptomyces')
  tlabs = tl$genus
#  tlabs[!tlabs %in% c(genus_root, genus1, genus2, genus3)] = ''
  tree$tip.label = tlabs
  
  tip.cols = rep('black', length(ids))
#  tip.cols[tl$genus == ] = 'seagreen'
  tip.cols[tl$genus %in% genus1] = 'lightpink2'
  tip.cols[tl$genus %in% genus2] = 'hotpink2'
  tip.cols[tl$genus %in% genus3] = 'deeppink3'
  
  scores = as.numeric(tree$node.label)
  node.bg = rep('white', tree$Nnode)
  node.bg[scores >= 0.9] = 'black'
  node.bg[scores >= 0.8 & scores < 0.9] = 'grey'
  tree$node.label = sprintf("%.02f", as.numeric(tree$node.label))
  
  fo = file.path(dirw, "19.pdf")
  pdf(fo, width = 11, height = 11)
  plot(tree, type = 'radial', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols, label.offset = 0.01, 
    no.margin = T, cex = 0.7)
#  tiplabels(pch = 22, frame = 'none', adj = 0.55, bg = tip.cols)
  nodelabels(pch = 22, bg = node.bg, cex = 0.5)
  add.scale.bar(x = 0.02, y = tree$Nnode*0.2 , lcol = 'black')
  
  legend(-0.3, 0.5, title = "# of Lantibiotics Characterized", legend = c("1 ", "3", "9"), fill = c("lightpink2", "hotpink2", "deeppink3"), xjust = 0.5, box.lwd = 0, cex=0.8)
  dev.off()


labs = sapply(strsplit(tl$taxa, split=";"), "[", 3)
tree$tip.label = labs
  fo = file.path(dirw, "19.lv3.pdf")
  pdf(fo, width = 11, height = 11)
  plot(tree, type = 'radial', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols, label.offset = 0.01, 
    no.margin = T, cex = 0.7)
  nodelabels(pch = 22, bg = node.bg, cex = 0.5)
  add.scale.bar(x = 0.02, y = tree$Nnode*0.2 , lcol = 'black')
  dev.off()
  
labs = sapply(strsplit(tl$taxa, split=";"), "[", 4)
tree$tip.label = labs
  fo = file.path(dirw, "19.lv4.pdf")
  pdf(fo, width = 11, height = 11)
  plot(tree, type = 'radial', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols, label.offset = 0.01, 
    no.margin = T, cex = 0.7)
  nodelabels(pch = 22, bg = node.bg, cex = 0.5)
  add.scale.bar(x = 0.02, y = tree$Nnode*0.2 , lcol = 'black')
  dev.off()
  