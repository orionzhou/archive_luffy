require(rtracklayer)
source("comp.plot.fun.R")

tname = "hm101"
qnames = c("hm004", "hm010", "hm018", "hm022", "hm034", "hm050", 
  "hm056", "hm058", "hm060", "hm095", "hm122", "hm129", 
  "hm185", "hm324", "hm340")
qnames = c("hm004", "hm010", "hm018", "hm034", "hm050", 
  "hm056", "hm058", "hm060", "hm095", "hm129", 
  "hm185", "hm324", "hm340")
cfgs = get_comp_cfg(tname, qnames)

dir = file.path(Sys.getenv("misc3"), 'compstat', 'figs')
fl = file.path(dir, 'loci.tbl')
tl = read.table(fl, header = T, sep = "\t", as.is = T)

i = 21
tls = tl[tl$i == i,]
gr =  GRanges(seqnames = tls$chr, ranges = IRanges(tls$beg, end = tls$end))
fn = sprintf("%s/fig%03d.pdf", dir, i)

source("comp.plot.fun.R")
dats = prep_plot_data(gr, cfgs, tname, qnames)


source("comp.plot.fun.R")
comp.plot(fn, dats, tname, qnames, 1000, subtitle = "pz")

