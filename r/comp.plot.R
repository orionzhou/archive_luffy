require(rtracklayer)
require(xlsx)
source("comp.plot.fun.R")

tname = "HM101"
qnames_all = c(
  "HM058", "HM125", "HM056", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340"
)
cfgs = get_comp_cfg(tname, qnames_all)

##### experimental
dir = file.path(Sys.getenv("misc3"), 'compstat', 'figs')
fl = file.path(dir, 'loci.xlsx')
tl = read.xlsx(fl, sheetIndex = 1, header = T)

for (i in 80:80) {
  #i = 1004
  tls = tl[tl$i == i,]
  
  qnames = sprintf("HM%03d", c(58, 129, 185, 4, 23, 340))
  
  gr =  GRanges(seqnames = tls$chr, ranges = IRanges(tls$beg, end = tls$end))
  fn = sprintf("%s/fig%03d.pdf", dir, i)

  dats = prep_plot_data(gr, cfgs, tname, qnames)
  comp.plot(fn, dats, tname, qnames, 1000, subtitle = "")
}

##### production
dir = file.path(Sys.getenv("misc3"), 'compstat', 'figs')
fl = file.path(dir, 'loci.xlsx')
tl = read.xlsx(fl, sheetIndex = 2, header = T)

for (i in 1:nrow(tl)) {
  tls = tl[i,]

  nums = strsplit(as.character(tls$orgs), split = ' ')[[1]]
  qnames = sprintf("HM%03d", as.numeric(nums))

  gr =  GRanges(seqnames = tls$chr, ranges = IRanges(tls$beg, end = tls$end))
  fn = sprintf("%s/r.%03d.pdf", dir, i)

  dats = prep_plot_data(gr, cfgs, tname, qnames)
  comp.plot(fn, dats, tname, qnames, 1000, subtitle = "")
}