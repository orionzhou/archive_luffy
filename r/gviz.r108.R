require(Gviz)
require(rtracklayer)
source("comp.fun.R")
source("gviz.fun.R")

dirw = file.path(Sys.getenv("misc2"), 'pb.stat')

org = "HM340.FN"
chr = 'scf005'
beg = 11094325
end = 12064779

org = "HM034"
chr = 'scf0000'
beg = 3235216
end = 4392933

cfg = get_genome_cfg(org)

	tgap = read.table(cfg$gap, header = F, sep = "\t", as.is = T)
	colnames(tgap) = c('chr','beg','end')
	tgap$beg = tgap$beg + 1
	tlen = read.table(cfg$size, header = F, sep = "\t", as.is = T)
	tlen = cbind(tlen, V3=1)
	colnames(tlen) = c('chr','end','beg')

axisTrack <- GenomeAxisTrack(cex = 1.0, exponent = 3)
ideoTrack <- build_ideogram_track(tgap, tlen, org)

gapTrack <- AnnotationTrack(genome = org, 
  chromosome = tgap$chr, start = tgap$beg, width = tgap$end-tgap$beg+1, 
  name = 'gap', showId = F, 
  fill = 'maroon', background.title = "midnightblue")

tg = read.table(cfg$gene, header = F, sep = "\t", as.is = T)
colnames(tg) = c("chr", "beg", "end", "srd", "id", "type", "cat")
tg = tg[tg$type == 'cds',]

dfg = data.frame(chromosome = tg$chr, start = tg$beg, end = tg$end,
  width = tg$end-tg$beg+1, strand = tg$srd, feature = tg$type, gene = tg$id, 
  exon = NA, transcript = tg$id, symbol = tg$id, stringsAsFactors = F)
gr.fill = rep('gray', nrow(dfg))
gr.fill[tg$cat %in% 'TE'] = 'orange'
gr.fill[tg$cat %in% 'Unknown'] = 'purple'
grTrack <- GeneRegionTrack(dfg, genome = org, shape = 'smallArrow',
  name = "genes", showId = F, just.group = 'below', stackHeight = 0.75,
  fill = gr.fill, fontsize = 9, max.height = 10, stacking = 'squish',
  cex.group = 0.8, cex.title = 1, background.title = 'midnightblue')


fn = sprintf("%s/%s_HM101/41_novseq/21.bed", Sys.getenv("misc3"), org)
tn = read.table(fn, header = F, sep = "\t", as.is = T)
colnames(tn) = c('chr','beg','end', 'type')
dfg = data.frame(chromosome = tn$chr, start = tn$beg, end = tn$end, stringsAsFactors = F)
nvTrack <- AnnotationTrack(dfg[dfg$chr==chr,], genome = org, name = 'Novel Segments', stacking = 'squish', background.title = 'tomato')


fo = sprintf("%s/41.%s.pdf", dirw, org)
pdf(file = fo, width = 8, height = 3, bg = 'transparent')
  plotTracks(
    list(axisTrack, gapTrack, nvTrack, grTrack),
    chromosome = chr, from = beg, to = end,
    min.height = 0, coverageHeight = 0.08, minCoverageHeight = 0,
    sizes = c(1, 0.4, 1, 1)
  )
dev.off()
