require(rtracklayer)
source("comp.plot.fun.R")

tname = "hm101"
qnames1 = c(
  "HM058", "HM125", "HM056", "HM129", "HM060", 
  "HM095", "HM185", "HM034", "HM004", "HM050", 
  "HM023", "HM010", "HM022", "HM324", "HM340"
)
qnames2 = c(
  "HM058", "HM129", 
  "HM185", "HM004", 
  "HM023", "HM340"
)
qnames = qnames2

cfgs = get_comp_cfg(tname, qnames)

dir = file.path(Sys.getenv("misc3"), 'compstat', 'figs')
fl = file.path(dir, 'loci.tbl')
tl = read.table(fl, header = T, sep = "\t", as.is = T)

for (i in 1003:1003) {
#i = 1004
tls = tl[tl$i == i,]
gr =  GRanges(seqnames = tls$chr, ranges = IRanges(tls$beg, end = tls$end))
fn = sprintf("%s/fig%03d.pdf", dir, i)

source("comp.plot.fun.R")
dats = prep_plot_data(gr, cfgs, tname, qnames)


source("comp.plot.fun.R")
comp.plot(fn, dats, tname, qnames, 1000, subtitle = "pz")

}