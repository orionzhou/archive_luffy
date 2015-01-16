require(plyr)
require(dplyr)
require(rtracklayer)
require(GenomicRanges)
require(ggplot2)
require(grid)
require(RColorBrewer)
source("comp.fun.R")

qnames = c(
  "HM058", "HM125", "HM056", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340"
)
tname = "HM101"
dir = file.path(Sys.getenv('misc3'), 'panseq')

##### refine mugsy result
fi = file.path(dir, '22.tbl')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
#ddply(ti, .(org), summarise, len = sum(end - beg + 1))

to = data.frame()
cid_max = max(ti$cid)
for (qname in qnames) {
#qname = "HM004"
idxs1 = which(ti$org == qname)
t1 = ti[idxs1,]
#t1 = t1[order(t1$chr, t1$beg),]

gr1 = with(t1, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
#grl1 = with(t1, makeGRangesListFromFeatureFragments(seqnames = chr, 
#  fragmentStarts = lapply(beg, rep), fragmentEnds = lapply(end, rep), 
#  strand = rep("*", nrow(t1))))

c1 = coverage(gr1)
l1 = mapply(function(x) {
  matrix(c(start(x), end(x), runValue(x)), nrow = 3, byrow = T)
}, c1)
lens = mapply(ncol, l1)
chrids = rep(names(l1), lens)
m1 = matrix(unlist(l1, use.names = F), ncol = 3, byrow = T)
colnames(m1) = c('beg', 'end', 'val')
t2 = data.frame(m1, stringsAsFactors = F)
t2 = cbind(chr = chrids, t2)
t3 = t2[t2$val > 1,]

idxs_rm = c()
for (i in 1:nrow(t3)) {
  chr = t3$chr[i]; beg = t3$beg[i]; end = t3$end[i]
  idxs = which(t1$org == qname & t1$chr == chr & t1$beg <= end & t1$end >= beg)
  stopifnot(length(idxs) > 1)
  lens = t1$end[idxs] - t1$beg[idxs] + 1
  idx_max = which(lens == max(lens))[1]
  idxs_rm = c(idxs_rm, idxs[-idx_max])
}
t2 = t1[-idxs_rm, ]
gr2 = with(t2, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
stopifnot(sum(width(gr2)) == sum(width(reduce(gr2))))

fb = sprintf("%s/%s_%s/41_novseq/21.bed", Sys.getenv("misc3"), qname, tname)
tb = read.table(fb, sep = "\t", header = F, as.is = T)
colnames(tb) = c('chr', 'beg', 'end', 'type')
tb = within(tb, {beg = beg + 1; chr = sprintf("%s_%d_%d", chr, beg, end);
  end = end - beg + 1; beg = 1})
grb = with(tb, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))

grr = setdiff(grb, gr2)
t2c = data.frame(cid = cid_max+1:length(grr), alen = width(grr), org = qname,
  chr = seqnames(grr), beg = start(grr), end = end(grr), srd = "+", 
  stringsAsFactors = F)
t4 = rbind(t2, t2c)
gr4 = with(t4, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
stopifnot(sum(width(gr4)) == sum(width(grb)))
cat(qname, ": ", sum(width(grr)), " overlap removed\n", sep = '')
to = rbind(to, t4)
}

fo = file.path(dir, "31.refined.tbl")
write.table(to, file = fo, sep = "\t", row.names = F, col.names = T, quote = F)

##### pan-16 AFS
fi = file.path(dir, '31.refined.tbl')
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ddply(ti, .(org), summarise, len = sum(end - beg + 1))

tu1 = ddply(ti, .(cid), summarise, n_org = length(org), 
  orgs = paste(sort(unique(as.character(org))), collapse = "_"),
  size = sum(end - beg + 1)
)

orgs = c(tname, qnames)
tu2 = ddply(tu1, .(n_org), summarise, size = sum(size))
dt1 = data.frame(n_org = tu2$n_org, size = tu2$size/tu2$n_org, org = 'mixed', stringsAsFactors = F)
dt1 = dt1[dt1$n_org != 1,]

tus = tu1[tu1$n_org == 1,]
dtt = ddply(tus, .(orgs), summarise, size = sum(size))
dt2 = data.frame(n_org = 1, size = dtt$size, org = dtt$orgs, stringsAsFactors = F)

dt = rbind(dt1, dt2)

cols = c(brewer.pal(11, 'Set2'), brewer.pal(4, 'Set1'), 'gray30')
labs = orgs

dt$org = factor(dt$org, levels = c(orgs, 'mixed'))
dt$n_org = factor(dt$n_org, levels = sort(as.numeric(unique(dt$n_org))))
p = ggplot(dt, aes(x = n_org, y = size/1000000, fill = org, order = plyr:::desc(org))) +
  geom_bar(stat = 'identity', position = "stack", width = 0.5) +
  scale_fill_manual(name = "Accession-Specific", breaks = labs, labels = labs, values = cols, guide = guide_legend(ncol = 1, byrow = F, label.position = "right", direction = "vertical", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  scale_x_discrete(name = '# Sharing Accession') +
  scale_y_continuous(name = 'Sequences (Mbp)', expand = c(0.02, 0)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(legend.position = "right", legend.key.size = unit(0.6, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line"), plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(axis.title.y = element_text(colour = 'black')) +
  theme(axis.text.x = element_text(size = 7, colour = "blue")) +
  theme(axis.text.y = element_text(size = 7, colour = "brown", angle = 0, hjust = 1))

fp = sprintf("%s/81.pan.afs.pdf", dir)
ggsave(p, filename = fp, width = 6, height = 5)

##### pan/core-genome size
fi = file.path(dir, '31.refined.tbl')
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

tp$rep = factor(tp$rep, levels = 1:max(tp$rep))
p = ggplot(tp) +
  geom_point(aes(x = n_org, y = pan/1000000), shape = 1, size = 1.1) +
  geom_point(aes(x = n_org, y = core/1000000), shape = 4, size = 1.1) +
#  geom_text(aes(x = n_org, y = 0, label = org), geom_params=list(size = 2.5, vjust = 0, angle = 30)) +
  stat_smooth(aes(x = n_org, y = pan/1000000, col = 'pan'), fill = 'azure4', size = 0.2) +
  stat_smooth(aes(x = n_org, y = core/1000000, col = 'core'), fill = 'azure4', size = 0.2) +
#  scale_shape(name = "", solid = FALSE, guide = F) +
  scale_color_manual(name = "", labels = c('Core-genome', 'Pan-genome'), values = c("dodgerblue", "firebrick1"), guide = guide_legend(label.position = "left", direction = "horizontal", title.theme = element_text(size = 8, angle = 0), label.theme = element_text(size = 8, angle = 0))) +
  scale_x_continuous(name = '# Genomes Sequenced') +
  scale_y_continuous(name = 'Sequences in Mbp', expand = c(0.02, 0), limits = c(0, 460)) + 
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.2, 'lines')) +
  theme(legend.position = "top", legend.key.size = unit(1, 'lines'), legend.background = element_rect(fill = 'white', size=0), legend.margin = unit(0, "line")) +
  theme(plot.margin = unit(c(0,1,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, color = "blue")) +
  theme(axis.text.y = element_text(size = 8, color = "grey", angle = 90, hjust  = 0.5))

fp = sprintf("%s/82.pan.size.pdf", dir)
ggsave(p, filename = fp, width = 6, height = 5)


# total NR segments, shared, accession-specific
fg = file.path(dir, '27.cluster.tbl')
tg = read.table(fg, header = T, sep = "\t", as.is = T)
tg = tg[tg$len >= 50,]
nrow(tg)
sum(tg$len)

tgs = tg[! tg$orgs %in% c("HM034", "HM056", "HM340.APECCA"),]
nrow(tgs)
sum(tgs$len)

sum(tg$len) - sum(tgs$len)

# sharing status of novel (coding) segments
fc = file.path(dir, '27.coord.tbl')
tc = read.table(fc, header = T, sep = "\t", as.is = T)
tc = cbind(tc, len = tc$end - tc$beg + 1)
sum(tc$len)

tc = tc[tc$len >= 50,]
ddply(tg, .(orgs), summarise, total_len = sum(len))
ddply(tc, .(orgs), summarise, total_len = sum(len))


fc034 = file.path(DIR_Data, 'genome', "HM034/51.bed/cds.bed")
fc056 = file.path(DIR_Data, 'genome', "HM056/51.bed/cds.bed")
fc340 = file.path(DIR_Data, 'genome', "HM340.APECCA/51.bed/cds.bed")
tc034 = read.table(fc034, sep = "\t", header = F, as.is = T)
tc056 = read.table(fc056, sep = "\t", header = F, as.is = T)
tc340 = read.table(fc340, sep = "\t", header = F, as.is = T)
colnames(tc034) = c("id", "beg", "end", 'gene', 'srd')
colnames(tc056) = c("id", "beg", "end", 'gene', 'srd')
colnames(tc340) = c("id", "beg", "end", 'gene', 'srd')
gr034 = GRanges(seqnames = tc034$id, ranges = 
  IRanges(tc034$beg+1, end = tc034$end))
gr056 = GRanges(seqnames = tc056$id, ranges = 
  IRanges(tc056$beg+1, end = tc056$end))
gr340 = GRanges(seqnames = tc340$id, ranges = 
  IRanges(tc340$beg+1, end = tc340$end))

org = "HM340"
orgs = "HM340 HM034 HM056"
tcs = tc[tc$org == org & tc$orgs == orgs,]
gr = GRanges(seqnames = tcs$id, ranges = IRanges(tcs$beg, end = tcs$end))

if(org == "HM034") {
  grc = gr034
} else if (org == "HM056") {
  grc = gr056
} else {
  grc = gr340
}
sum(width(gr))
sum(width(intersect(gr, grc)))

# plot novel segments length distribution
tmp1 = table(tc$len)
tp.1 = data.frame(len=as.numeric(names(tmp1)), cnt=c(tmp1))
tp.2 = cbind(tp.1, sum=tp.1$len * tp.1$cnt)
tp.3 = tp.2[order(tp.2$len, decreasing=T),]
tp = cbind(tp.3, cumsum=cumsum(tp.3$sum))
for (x in c(1,50,60,100,500,1000)) {
  tt = tc[tc$len >= x,]
  cat(x, nrow(tt), sum(tt$len), "\n", sep = "\t")
}
#segments(0, y, x, y, col='blue')
#segments(x, y, x, 0, col='blue')
