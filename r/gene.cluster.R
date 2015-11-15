require(grid)
require(plyr)
require(dplyr)
require(ggplot2)
source('Location.R')
source("comp.fun.R")

dirw = file.path(Sys.getenv("misc2"), "gene.cluster")
dir.create(dirw)

org = "HM101"

f_gen = file.path(Sys.getenv("genome"), org, "51.gtb")
tg = read.table(f_gen, header = T, sep = "\t", as.is = T)[,c(1,3:5,16:17)]

### process mcl output
fi = sprintf("%s/04.mcl/%s.mcl", dirw, org)
ti = read.table(fi, sep = "-", header = F, as.is = T)
ti = cbind(grp = 1:nrow(ti), ti)
x = apply(ti, 1, kkk <- function(rw) { 
  data.frame(grp = as.numeric(rw[1]), id = unlist(strsplit(rw[2], "\t")), stringsAsFactors = F)
})
require(data.table)
tc = rbindlist(x)
tc = merge(tg, tc, by = 'id', all.x = T)
y = which(is.na(tc$grp))
tc$grp[y] = seq(from = max(tc$grp, na.rm = T)+1, by = 1, length.out = length(y))

tt = tc[tc$cat2 != 'TE',]
tt = tt[order(tt$chr, tt$beg, tt$end),]
tt = cbind(idx = 1:nrow(tt), tt)
to = merge(tc[,c('id','grp')], tt[,c('idx','id')], by = 'id', all.x = T)

fo = sprintf("%s/08.grp/%s.tbl", dirw, org)
write.table(to, fo, col.names = T, row.names = F, sep = "\t", quote = F)


### identify tandem clusters
fi = sprintf("%s/08.grp/%s.tbl", dirw, org)
ti = read.table(fi, sep = "\t", header = T, as.is = T)
ti = ti[!is.na(ti$idx),]
ti = ti[order(ti$idx),]

aids = c(); bids = c()
wsize = 2
for (i in 1:nrow(ti)) {
  for (j in (i+1):(i+wsize)) {
    if(j > nrow(ti)) next
    if(ti$grp[i] == ti$grp[j]) {
      aids = c(aids, ti$id[i])
      bids = c(bids, ti$id[j])
    }
  }
}

tt = data.frame(aid = aids, bid = bids, stringsAsFactors = F)
write.table(tt, 'tmp.tbl', col.names = T, row.names = F, sep = "\t", quote = F)
system("graph.connComp.py tmp.tbl tmpo.tbl")
tc = read.table("tmpo.tbl", sep = "\t", header = T, as.is = T)
system(sprintf("rm %s %s", "tmp.tbl", "tmpo.tbl"))

to = merge(tg, tc, by = 'id', all.x = T)
fo = sprintf("%s/11.tandem/%s.tbl", dirw, org)
write.table(to, fo, col.names = T, row.names = F, sep = "\t", quote = F, na = '')

### cluster distribution for different families
fi = sprintf("%s/11.tandem/%s.tbl", dirw, org)
ti = read.table(fi, header = T, sep = "\t", as.is = T)
ti = ti[ti$cat2 != 'TE',]
idx_single = which(is.na(ti$clu))
ti$clu[idx_single] = seq(max(ti$clu, na.rm = T)+1, by = 1, length.out = length(idx_single))

gb = group_by(ti, clu)
tc = summarise(gb, csize = n())
ti = merge(ti, tc, by = 'clu')

brks = c(seq(0.5, 10.5, by = 1), 15.5, Inf)
labs = c(1:10, '11-15', '16+')
labs = factor(labs, levels = labs)

to = ti
to$cat2[to$cat2 %in% c("CC-NBS-LRR", "TIR-NBS-LRR")] = "NBS-LRR"
fams = c("NBS-LRR", "F-box", "LRR-RLK", "NCR", "Unknown", "CRP0000-1030", "CRP1600-6250")
to$cat2[! to$cat2 %in% fams] = 'Pfam-Miscellaneous'

do = data.frame()
for (fam in c('NCR', 'NBS-LRR', 'LRR-RLK', "F-box", 'Pfam-Miscellaneous')) {
  tm = to[to$cat2 == fam,]
  itvs = table(cut(tm$csize, breaks = brks))
  ds = data.frame(fam = fam, itv = labs, n = as.numeric(itvs), prop = as.numeric(itvs) / sum(itvs), stringsAsFactors = T)
  do = rbind(do, ds)
}

x = do[do$itv == '1',]
fams = as.character(x$fam[order(x$prop, decreasing = T)])
do$fam = factor(do$fam, levels = fams)

p1 = ggplot(do) + 
  geom_bar(aes(x = itv, y = prop), stat = 'identity', geom_params = list(width = 0.8)) + 
  scale_x_discrete(name = 'Tandem array size') +
  scale_y_continuous(name = 'Proportion in family') +
  facet_wrap(~ fam, scales = 'free', nrow = 2) +  
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, angle = 60, colour = "blue", hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 0.5))

fo = sprintf("%s/12.tandem.plot/%s.pdf", dirw, org)
ggsave(p1, filename = fo, width = 7, height = 5)


##### tandem array stats in  different genomes
do = data.frame()

for (org in orgs) {
#org = "HM101"
fi = sprintf("%s/12_tandem/%s.tbl", dirw, org)
ti = read.table(fi, header = F, sep = "\t", as.is = T)
colnames(ti) = c('clu', 'id')

x = ddply(ti, .(clu), summarise, size = length(id))
itvs = table(cut(x$size, breaks = c(seq(1.5, 10.5, by = 1), 15.5, 20.5, Inf)))
labs = c(2:10, '11-15', '16-20', '21+')
labs = factor(labs, levels = labs)

ds = data.frame(org = org, itv = labs, n = as.numeric(itvs), stringsAsFactors = T)
do = rbind(do, ds)
}
do$org = factor(do$org, levels = orgs)

p1 = ggplot(do) + 
  geom_bar(aes(x = itv, y = n, fill = org),
    stat = 'identity', position = 'dodge', geom_params = list(width = 0.8)) + 
  scale_fill_manual(name = "", values = c("#E41A1C", "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FDBF6F","#FF7F00")) +
  scale_x_discrete(name = 'Size of tandem gene arrays') +
  scale_y_continuous(name = 'Number of arrays') +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines'), axis.ticks.margin = unit(0.4, 'lines')) +
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank(), legend.background = element_rect(fill = 'white', colour = 'black', size = 0.3), legend.key = element_rect(fill = NA, colour = NA, size = 0), legend.key.size = unit(1, 'lines'), legend.margin = unit(0, "lines"), legend.title = element_text(size = 9, angle = 0), legend.text = element_text(size = 9, angle = 0)) +
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(axis.title.x = element_text(size = 9, angle = 0)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "blue")) +
  theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 90, hjust = 1))

fo = sprintf("%s/51.pdf", dirw)
ggsave(p1, filename = fo, width = 7, height = 5)

for (org in orgs) {

f_gen = file.path(Sys.getenv("genome"), org, "51.gtb")
tg = read.table(f_gen, header = T, sep = "\t", as.is = T)[,c(1,3:5,15:17)]

fi = sprintf("%s/12_tandem/%s.tbl", dirw, org)
ti = read.table(fi, header = F, sep = "\t", as.is = T)
colnames(ti) = c('clu', 'id')

pct = nrow(ti)/nrow(tg)
cat(sprintf("%10s: %5d genes, %d clusters | %d (%.03f) in tandem arrays\n", org, nrow(tg), length(unique(ti$clu)), nrow(ti), pct))
}
