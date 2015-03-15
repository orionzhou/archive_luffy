require(rtracklayer)
require(xlsx)
source("comp.fun.R")
source("comp.plot.fun.R")

dirw = file.path(Sys.getenv("misc3"), 'comp.stat', 'figs')

##### experimental
source("comp.plot.fun.R")
fl = file.path(dirw, 'loci.xlsx')
tl = read.xlsx(fl, sheetIndex = 1, header = T)

tracks = c('qgene', 'qaxis', 'qgap', 'link', 'tgap', 'taxis', 'tgene')
#for (i in 80:80) {
  i = 16
  tls = tl[tl$i == i,]
  
  qnames = sprintf("HM%03d", c(58, 34, 129, 185, 4, 23, 340))
#  qnames = qnames_all
  
  gro =  GRanges(seqnames = tls$chr, ranges = IRanges(tls$beg, end = tls$end))
  fn = sprintf("%s/fig%03d.pdf", dirw, i)

dats = prep_plot_data(gro, cfgs, tname, qnames, tracks)
res = comp.plot(dats, tname, qnames, tracks, draw.title = F)

pdf(file = fn, width = 7, height = res$ht/72, bg = 'transparent')
grid.newpage()
grid.draw(res$grobs)
dev.off()
#}

##### production
fl = file.path(dirw, 'loci.xlsx')
tl = read.xlsx(fl, sheetIndex = 2, header = T)

for (i in 1:nrow(tl)) {
  tls = tl[i,]

  nums = strsplit(as.character(tls$orgs), split = ' ')[[1]]
  qnames = sprintf("HM%03d", as.numeric(nums))

  gr =  GRanges(seqnames = tls$chr, ranges = IRanges(tls$beg, end = tls$end))
  fn = sprintf("%s/r.%03d.pdf", dirw, i)

  dats = prep_plot_data(gr, ccfgs, tname, qnames)
  comp.plot(fn, dats, tname, qnames, 1000, subtitle = "")
}


### tandem duplication illustration
chr = "chr8"
beg = 3370714
end = 3416737
gro =  GRanges(seqnames = chr, ranges = IRanges(beg, end = end))
source("comp.plot.fun.R")

qnames = sprintf("HM%03d", c(56, 10, 129, 34))
qnames = c("HM034")
tracks = c('tgene', 'trnaseq', 'taxis', 'tgap', 'link', 'qgap', 'qaxis', 'qgene', 'qpacbio', 'qrnaseq')
  
dats = prep_plot_data(gro, cfgs, tname, qnames, tracks)
res = comp.plot(dats, tname, qnames, tracks, draw.title = T)

fn = sprintf("%s/illus_tandup.pdf", dirw)
pdf(file = fn, width = 7, height = res$ht/72, bg = 'transparent')
grid.newpage()
grid.draw(res$grobs)
dev.off()

### SNP calling inconsistency illustration
chr = "chr5"
beg = 1555000
end = 1605000
gro =  GRanges(seqnames = chr, ranges = IRanges(beg, end = end))
source("comp.plot.fun.R")

qnames = c("HM125")
tracks = c('qgene', 'qaxis', 'qgap', 'link', 'tgap', 'taxis', 'tgene', 'tmapp', 'mcov', 'msnp', 'tsnp')
  
dats = prep_plot_data(gro, cfgs, tname, qnames, tracks)
res = comp.plot(dats, tname, qnames, tracks, draw.title = F)

fn = sprintf("%s/illus_snpcall.pdf", dirw)
pdf(file = fn, width = 7, height = res$ht/72, bg = 'transparent')
grid.newpage()
grid.draw(res$grobs)

xs = c(88, 250, 407)
wds = c(155, 63, 35)
ys = rep(3, 3)
hts = rep(103, 3)
grid.rect(x = unit(xs, 'points'), width = unit(wds, 'points'), 
  y = unit(ys, 'points'), height = unit(hts, 'points'),
  just = c('left', 'bottom'),
  gp = gpar(lwd = 1, col = 'dodgerblue3', fill = NA)
)
labs = LETTERS[1:length(xs)]
grid.text(labs, x = unit(xs, 'points'), y = unit(ys+hts, 'points'), 
  just = c('left', 'top'), 
  gp = gpar(col = "dodgerblue3", fontface = 2, fontsize = 13)
)
dev.off()

##### SV illustration
### look for del/ins/inv/tlc
fb = '/home/youngn/zhoup/Data/misc3/HM034_HM101/31_sv/05.stb'
tb = read.table(fb, header = T, sep = "\t", as.is = T)
fx = '/home/youngn/zhoup/Data/misc3/HM034_HM101/31_sv/05.stx'
tx = read.table(fx, header = T, sep = "\t", as.is = T)

# ins/del
td = tx[tx$type == "DEL",]
td[td$tend - td$tbeg > 5000,][21:40,]
tn = tx[tx$type == "INS",]
tn[tn$qend - tn$qbeg > 3000,][1:20,]

# tlc
tc = tx[tx$type %in% c('TLC:INS', 'TLC:DEL'),]
tb = tb[,c(1:4,8:10)]
to = merge(tc, tb, by = 'id')
to = to[order(to$tchr.x, to$tbeg.x),]
to[to$tend.x-to$tbeg.x>10000,][41:60,]

# look for inv
idxs = c()
for (i in 1:(nrow(to)-1)) {
  if(to$id[i] == to$id[i+1] & to$tbeg.x[i] == to$tbeg.x[i+1]) {
    idxs = c(idxs, i)
  }
}

fl = file.path(dirw, 'loci.xlsx')
tl = read.xlsx(fl, sheetIndex = 3, header = T)

idxs = sort(unique(tl$idx))
tracks = c('qgene', 'qaxis', 'qgap', 'link', 'tgap', 'taxis', 'tgene')
ress = list()
for (idx in idxs) {
  tls = tl[tl$idx == idx,]

  nums = strsplit(as.character(tls$org), split = ' ')[[1]]
  qnames = sprintf("HM%03d", as.numeric(nums))

  gro =  GRanges(seqnames = tls$chr, ranges = IRanges(tls$beg, end = tls$end))

  dats = prep_plot_data(gro, cfgs, tname, qnames, tracks)
  res = comp.plot(dats, tname, qnames, tracks, draw.title = F)
  ress[[idx]] = res
}
ht = ress[[1]]$ht
htm = 15

nplot = length(idxs)
fn = sprintf("%s/illus_sv.pdf", dirw)
pdf(file = fn, width = 7, height = (ht+htm)*nplot/72, bg = 'transparent')
  
numrow = nplot*2; numcol = 3
grid.newpage()
  
pushViewport(viewport(layout = grid.layout(numrow, numcol, 
  heights = rep(c(htm, ht), nplot),
  widths = unit.c(unit(htm, 'points'), unit(1, 'npc') - unit(htm*2, 'points'), unit(htm, 'points')))
))

for (idx in idxs) {
  pushViewport(viewport(layout.pos.row = idx*2, layout.pos.col = 2))
  grid.draw(ress[[idx]]$grobs)
  popViewport()
}
labs = LETTERS[1:nplot]
grid.text(labs, x = unit(2, 'points'), 
  y = unit(c(nplot:1) * (ht+htm), 'points'), 
  just = c('left', 'top'), 
  gp = gpar(col = "black", fontface = 2, fontsize = 16)
)
dev.off()
