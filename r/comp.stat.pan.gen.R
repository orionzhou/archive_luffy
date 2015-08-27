require(plyr)
require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(RColorBrewer)
require(dplyr)
source("comp.fun.R")

dirw = file.path(Sys.getenv('misc3'), 'comp.panseq')
diro = file.path(Sys.getenv("misc3"), "comp.stat")

##### novel segments AFS
fi = file.path(dirw, '31.refined.tbl')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ddply(ti, .(org), summarise, len = sum(end - beg + 1))

gb = group_by(ti, cid)
tu1 = summarise(gb, n_org = n(), 
  orgs = paste(sort(unique(as.character(org))), collapse = "_"),
  size = sum(end - beg + 1))

orgs = c(tname, qnames)
tu2 = ddply(tu1, .(n_org), summarise, size = sum(size))
dt1 = data.frame(n_org = tu2$n_org, size = tu2$size/tu2$n_org, org = 'mixed', stringsAsFactors = F)
dt1 = dt1[dt1$n_org != 1,]

tus = tu1[tu1$n_org == 1,]
dtt = ddply(tus, .(orgs), summarise, size = sum(size))
dt2 = data.frame(n_org = 1, size = dtt$size, org = dtt$orgs, stringsAsFactors = F)
dt2$org = factor(dt2$org, levels = orgs)
dt2 = dt2[order(dt2$org),]

to = rbind(dt2, dt1)
mxs = seq(20, length.out = 14, by = 5)
to = cbind(to, x = c(1:15, mxs), wid = c(rep(1, 15), rep(2, 14)))

cols = c(brewer.pal(11, 'Set3'), brewer.pal(4, 'Set1'), 'gray30')
labs = orgs

#to$org = factor(to$org, levels = c(orgs, 'mixed'))
#to$n_org = factor(to$n_org, levels = sort(as.numeric(unique(to$n_org))))
p1 = ggplot(to, aes(x = x, y = size/1000000, fill = org, width = wid)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(name = "Accession-Specific", breaks = labs, labels = labs, values = cols) +
  scale_x_continuous(name = '# Sharing Accession', breaks = c(8, mxs), labels = c('1(accession-specific)', 2:15), expand = c(0, 0), limits = c(-2, max(to$x) + 3)) +
  scale_y_continuous(name = 'Sequences (Mbp)', expand = c(0, 0), limits = c(0, 8.5)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(legend.position = c(0.7, 0.7), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(0.5, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_text(size = 8, angle = 0), legend.text = element_text(size = 8, angle = 0)) +
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))

#fp = sprintf("%s/61.pan.genome.afs.pdf", diro)
#ggsave(p, filename = fp, width = 5, height = 5)

##### pan/core-genome size curve
fi = file.path(dirw, '31.refined.tbl')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti = within(ti, {len = end - beg + 1})
ti = ti[ti$len > 20,] ##### filter very short segments

tcfg = get_genome_cfg(tname)
tt1 = read.table(tcfg$size, sep = "\t", as.is = T)
grt1 = with(tt1, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt2 = read.table(tcfg$gap, sep = "\t", as.is = T)
grt2 = with(tt2, GRanges(seqnames = V1, ranges = IRanges(V2 + 1, end = V3)))
grt = setdiff(grt1, grt2)

grgs = list()
for (qname in qnames) {
  fgax = sprintf("%s/%s_%s/23_blat/31.5/gax", Sys.getenv("misc3"), qname, tname)
  gax = read.table(fgax, header = F, sep = "\t", as.is = T)
  grg = with(gax, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
  grgs[[qname]] = grg
}

reps = 1:8
n_orgs = 1:16
tp = data.frame(rep = rep(reps, each = length(n_orgs)), 
  n_org = rep(n_orgs, length(rep)), core = NA, pan = NA)
  
for (rep in 1:8) {
  grc = grt
  grp = grt
  core = sum(width(grc))
  pan = sum(width(grc))
  tp$core[tp$rep == rep && tp$n_org == 1] = core
  tp$pan[tp$rep == rep && tp$n_org == 1] = pan
  
  set.seed(rep*100)
  qnames.rep = sample(qnames)
  for (i in 1:length(qnames.rep)) {
    qname = qnames.rep[i]
    org_str = paste(c(tname, qnames.rep[1:i]), collapse = "+")
    
    grc = intersect(grc, grgs[[qname]])
    core = sum(width(grc))
    tp$core[tp$rep == rep & tp$n_org == i+1] = core
    
    tis = ti[ti$org %in% qnames.rep[1:i], ]
#    tu1 = ddply(tis, .(cid), summarise, n_org = length(org), size = sum(len), .parallel = T)
    tis_df = group_by(tis, cid)
    tu1 = summarise(tis_df, n_org = length(org), size = sum(len))
#    system.time(summarise(tis_df, n_org = length(org), size = sum(len)))
    tu2 = ddply(tu1, .(n_org), summarise, size = sum(size))
    tu3 = cbind(tu2, persize = tu2$size / tu2$n_org)
    pan = sum(width(grp)) + sum(tu3$persize)
    tp$pan[tp$rep == rep & tp$n_org == i+1] = pan
    
    cat(rep, qname, core, pan, "\n")
  }
}
fo = file.path(diro, "63.pan.genome.size.tbl")
write.table(tp, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

##plot
fi = file.path(diro, "63.pan.genome.size.tbl")
tp = read.table(fi, header = T, sep = "\t", as.is = T)

tp$rep = factor(tp$rep, levels = 1:max(tp$rep))
p2 = ggplot(tp) +
  geom_point(aes(x = n_org, y = pan/1000000), shape = 1, size = 1) +
  geom_point(aes(x = n_org, y = core/1000000), shape = 4, size = 1) +
#  geom_text(aes(x = n_org, y = 0, label = org), geom_params=list(size = 2.5, vjust = 0, angle = 30)) +
  stat_smooth(aes(x = n_org, y = pan/1000000, col = 'a'), fill = 'azure4', size = 0.3, se = F) +
  stat_smooth(aes(x = n_org, y = core/1000000, col = 'b'), fill = 'azure4', size = 0.3, se = F) +
#  scale_shape(name = "", solid = FALSE, guide = F) +
  scale_color_manual(name = "", labels = c('Pan-genome', 'Core-genome'), values = c("firebrick1", "dodgerblue")) +
  scale_x_continuous(name = '# Genomes Sequenced') +
  scale_y_continuous(name = 'Genome size (Mbp)', expand = c(0, 0), limits = c(0, 480)) + 
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.2, 'lines')) +
  theme(legend.position = c(0.15, 0.2), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA), legend.key.size = unit(1, 'lines'), legend.margin = unit(0, "line"), legend.title = element_blank(), legend.text = element_text(size = 8, angle = 0)) +
  theme(plot.margin = unit(c(0,1,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, color = "blue")) +
  theme(axis.text.y = element_text(size = 8, color = "grey", angle = 90, hjust  = 0.5))

#fp = sprintf("%s/63.pan.genome.size.pdf", dir)
#ggsave(p, filename = fp, width = 5, height = 4)

fo = sprintf("%s/63.pan.genome.pdf", diro)
pdf(file = fo, width = 10, height = 5, bg = 'transparent')
numrow = 1; numcol = 2
grid.newpage()
pushViewport(viewport(layout = grid.layout(numrow, numcol)))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

dco = data.frame(x = rep(1:numrow, each = numcol), y = rep(1:numcol, numrow), lab = LETTERS[1:(numrow*numcol)])
for (i in 1:nrow(dco)) {
  x = dco$x[i]; y = dco$y[i]; lab = dco$lab[i]
  grid.text(lab, x = 0, y = unit(1, 'npc'), just = c('left', 'top'), gp = gpar(col = "black", fontface = 2, fontsize = 20),
    vp = viewport(layout.pos.row = x, layout.pos.col = y))
}
dev.off()


