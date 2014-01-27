require(ape)
plot_mt_tree_2 <- function(fi, fo, ann, opt) {
  tree = read.tree(fi)

  label_id = tree$tip.label
  idx_HM101 = which(label_id == 'HM101')
  idx_acc26 = which(label_id %in% get_mt_ids('acc26'))

  tip.color = rep('black', length(tree$tip.label))
  tip.color[idx_acc26] = 'royalblue'
  tip.color[idx_HM101] = 'red'
  
  df1 = data.frame(id=tree$tip.label)
  df2 = merge(df1, ann, by="id", all.x=TRUE)
  label_origin = paste(df2$country, df2$category, sep=" | ")
  tree$tip.label = label_origin

  scores = tree$node.label
  node.labels.bg = rep('white', tree$Nnode)
  node.labels.bg[scores >= 0.95] = 'black'
  node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

  if(opt == 'acc26') {
    cex1 = 1
    cex2 = 1.2
    label.offset = 0.07
    xmax = 0.5
  } else if (opt == 'acc31') {
    cex1 = 1
    cex2 = 1.2
    label.offset = 0.07
    xmax = 0.75
  } else {
    stop(cat("unknown opt: ", opt, "\n", sep=""))
  }

  png(filename=fo, width=600, height=600, units='px')
  plot(tree, show.node.label=FALSE, show.tip.label=TRUE, tip.color=tip.color, font=3, no.margin=TRUE, 
    cex=cex1, label.offset=label.offset, x.lim=xmax)
  tiplabels(label_id, adj=0, bg=NA, font=2, frame='none', col=tip.color, cex=cex2)
  nodelabels(pch = 21, bg=node.labels.bg)
  add.scale.bar()
  dev.off()
}
plot_mt_tree_3 <- function(fi, fo, ann) {
  tree = read.tree(fi)

  group1 = c("HM101", "HM056", "HM058", "HM117", "HM125")
  group2 = c("HM340", "HM324", "HM018", "HM022-I", "HM017-I")
  group3 = c("HM034", "HM129", "HM060", "HM095", "HM185")
  
  labels = tree$tip.label
  tip.color = rep('black', length(tree$tip.label))
  tip.color[which(labels %in% group1)] = 'red'
  tip.color[which(labels %in% group2)] = 'forestgreen'
  tip.color[which(labels %in% group3)] = 'dodgerblue'

  df1 = data.frame(idx = 1 : length(labels), id = labels)
  df2 = merge(df1, ann, by = "id", all.x = T)
  df3 = df2[order(df2$idx), ]
  labelsn = paste(df3$id, df3$country, df3$category, sep = " | ")
  tree$tip.label = labelsn

  scores = as.numeric(tree$node.label)
  if(mean(scores, na.rm = TRUE) > 1) { scores = scores / 1000 }
  node.labels.bg = rep('white', tree$Nnode)
  node.labels.bg[scores >= 0.95] = 'black'
  node.labels.bg[scores >= 0.8 & scores < 0.95] = 'gray'

  png(filename=fo, width=2000, height=5000, units='px')
  plot(tree, show.node.label = F, show.tip.label = T, tip.color = tip.color, 
    font = 2, no.margin = T, cex = 1.2)
  nodelabels(pch = 22, bg = node.labels.bg)
  add.scale.bar(x = 0.05, y = 20, lcol = 'black')
  dev.off()
}

f_ann = file.path(DIR_Data, "misc3/hapmap_mt40/31_phylogeny/mt_label.tbl")
ann = read.table(f_ann, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

opt = 'all'
reg = "chr5"
dir = file.path(DIR_Data, "misc3/hapmap_mt40/31_phylogeny", opt)
fi = sprintf("%s/22_phyml/%s.nwk", dir, reg)
fo = sprintf("%s/22_phyml/%s.png", dir, reg)
#fi = sprintf("%s/21_phynj/%s.phb", dir, reg)
#fo = sprintf("%s/21_phynj/%s.png", dir, reg)
plot_mt_tree_3(fi, fo, ann)


opt = 'acc31'
dir = file.path(DIR_Data, "repo/mt_35/31_phylogeny", opt)
for(i in 1:8) {
  fi = sprintf("%s/22_phyml/chr%d.nwk", dir, i)
  fo = sprintf("%s/22_phyml/chr%d.png", dir, i)
  plot_mt_tree(fi, fo, ann)
}
fi = file.path(dir, "19/04.nwk")
fo = file.path(dir, "19/04.png")
plot_mt_tree(fi, fo)

opt = 'acc56'
dir = file.path(DIR_Data, "repo/mt_35/31_phylogeny", opt)
for(i in 1:8) {
  fi = sprintf("%s/22_phyml/chr%d.nwk", dir, i)
  fo = sprintf("%s/22_phyml/chr%d.png", dir, i)
  plot_mt_tree(fi, fo, ann)
}


#compare sv phylogeny with chr5 phylogeny
  fi_chr5 = file.path(DIR_Data, "repo/mt_35/31_phylogeny", "acc26", "21_phynj/chr5.phb")
  fo_chr5 = file.path(DIR_Data, "repo/mt_35/31_phylogeny", "acc26", "21_phynj/chr5.png")
  tree = read.tree(fi_chr5)

  tip.text = tree$tip.label
  tip.bg.chr5 = rainbow(length(tree$tip.label))
  names(tip.bg.chr5) = tip.text

  scores = as.numeric(tree$node.label)
  if(mean(scores, na.rm=TRUE) > 1) { scores = scores / 1000 }
  node.bg = rep('white', tree$Nnode)
  node.bg[scores >= 0.95] = 'black'
  node.bg[scores >= 0.8 & scores < 0.95] = 'gray'

  png(filename=fo_chr5, width=500, height=500, units='px')
  plot(tree, show.node.label=FALSE, show.tip.label=FALSE, font=4, no.margin=TRUE, cex=1.2)
  tiplabels(text=tip.text, bg=tip.bg.chr5)
  nodelabels(pch=22, bg=node.bg)
  add.scale.bar()
  dev.off()

  fi_sv = file.path(DIR_Repo, "mt_35/40_sv/41_shared/73_phylogeny", "04.phb")
  fo_sv = file.path(DIR_Repo, "mt_35/40_sv/41_shared/73_phylogeny", "04.png")
  tree = read.tree(fi_sv)
  tree = rotate(tree, c("HM015", "HM003"))
#  tree = rotate(tree, c("HM101", "HM008"))
#  tree = rotate(tree, c("HM101", "HM007"))
  tree = rotate(tree, c("HM020", "HM021"))
  tree = rotate(tree, c("HM021", "HM028"))

  tip.text = tree$tip.label
  tip.bg.sv = c(tip.bg.chr5[tip.text])

  scores = as.numeric(tree$node.label)
  if(mean(scores, na.rm=TRUE) > 1) { scores = scores / 1000 }
  node.bg = rep('white', tree$Nnode)
  node.bg[scores >= 0.95] = 'black'
  node.bg[scores >= 0.8 & scores < 0.95] = 'gray'

  png(filename=fo_sv, width=500, height=500, units='px')
  plot(tree, show.node.label=FALSE, show.tip.label=FALSE, font=4, no.margin=TRUE, cex=1.2)
  tiplabels(text=tip.text, bg=tip.bg.sv)
  nodelabels(pch=22, bg=node.bg)
  add.scale.bar()
  dev.off()

  

