require(Gviz)
require(rtracklayer)
source("comp.fun.R")

tname = "hm101"
qname = "hm056"
rname = "hm340"

t = read_genome_stat(tname)
q = read_genome_stat(qname)
r = read_genome_stat(rname)

vq = read_var_stat(qname)
vr = read_var_stat(rname)

cq = read_comp_stat(qname, tname)
cr = read_comp_stat(rname, tname)

# build target-tracks
axisTrack <- GenomeAxisTrack(cex=1.0, exponent=3)
ideoTrack <- build_ideogram_track(t$gap, t$len, tname)

gapTrack <- AnnotationTrack(genome = tname, 
  chromosome = t$gap$id, start = t$gap$beg, width = t$gap$len, 
  name = 'gap', showId = F, 
  fill = 'maroon', background.title = "midnightblue")

f_mapp = "/home/youngn/zhoup/Data/genome/pan3/18_stat_k60/15_mapp.bw"
mappTrack <- DataTrack(range = f_mapp, genome = tname, type = 'h', 
  showAxis = F, name = 'mapp', background.title = 'midnightblue')

grTrack <- GeneRegionTrack(t$gene, genome = tname, 
  name = "genes", showId = T, just.group = 'below', stackHeight = 0.5,
  cex.group = 0.8, background.title = 'midnightblue')

# plot
fg = file.path(t$dir, "51.gtb")
tg = read.table(fg, header=T, sep="\t", quote="", as.is=T)[,c(1,3:6,15:17)]
tgs = tg[grepl('CRP', tg$cat3) | grepl('NBS', tg$cat3),]
tgs = tgs[tgs$chr %in% sprintf("chr%d", 1:8),]
tgs = tgs[order(tgs$id),]

for (i in 1:100) {
  chr = tgs$chr[i]
  beg = tgs$beg[i]
  end = tgs$end[i]
  fo = sprintf("%s/figs/genes/%s.png", cq$dir, tgs$id[i])

  qvar <- build_var_tracks(vq, chr, beg, end, qname, tname)
  rvar <- build_var_tracks(vr, chr, beg, end, rname, tname)
  qcomp <- build_comp_tracks(cq, chr, beg, end, qname, tname)
  rcomp <- build_comp_tracks(cr, chr, beg, end, rname, tname)

  CairoPNG(filename = fo, width = 900, height = 1200)
  plotTracks(
    list(ideoTrack, axisTrack, gapTrack, mappTrack, grTrack, 
      qcomp$compTrack, qcomp$snpTrack, qcomp$siTrack, qcomp$liTrack, 
      qvar$snpTrack, qvar$indelTrack, qvar$covTrack,
      rcomp$compTrack, rcomp$snpTrack, rcomp$siTrack, rcomp$liTrack,
      rvar$snpTrack, rvar$indelTrack, rvar$covTrack),
    chromosome = chr, from = beg, to = end, 
    extend.left = (end - beg) / 20, extend.right = (end - beg) / 20, 
    sizes = c(1, 1, 1, 2, 2,  
      2, 2, 2, 2, 
      1, 1, 1,  
      2, 2, 2, 2, 
      1, 1, 1))
  dev.off()
}

# comparison plot using a region
f_reg = '/home/youngn/zhoup/Data/genome/HM101/81_regions.tbl'
reg = read.table(f_reg, header=T, sep="\t", as.is=T)
i = 5
chr=reg$chr[i]; beg=reg$beg[i]; end=beg+30000


# dotplot (outdated)
xb = c()
xe = c()
yb = c()
ye = c()
for (i in 1:nrow(tbs)) {
  xb = c(xb, tbs$hBeg[i])
  xe = c(xe, tbs$hEnd[i])
  if(tbs$qSrd[i] == '-') {
    yb = c(yb, tbs$qEnd[i])
    ye = c(ye, tbs$qBeg[i])
  } else {
    yb = c(yb, tbs$qBeg[i])
    ye = c(ye, tbs$qEnd[i])
  }
}
tt = cbind(tbs, xb = xb, xe = xe, yb = yb, ye = ye)
p <- ggplot(tt) +
  geom_segment(mapping = aes(x = xb, xend = xe, y = yb, yend = ye), 
    size = 0.6) +
  layer(data = dat2.coord, geom = 'rect', 
    mapping = aes(xmin = hBeg.r, xmax = hEnd.r, ymin = -3000000, ymax = 0, 
      fill = hId), geom_params = list(size = 0)) +
  layer(data = dat2.coord, geom = 'text', 
    mapping = aes(x = (hBeg.r + hEnd.r) / 2, y = -4000000, label = hId), 
    geom_params = list(hjust = 0.5, vjust = 1, angle = 0, size = 3)) +
  facet_grid(qId ~ hId, scales = "free") +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(name = 'HM101 (Mt4.0)') +
  scale_y_continuous(name = dat1.name) +
  theme_bw()
  theme(legend.position = "none")
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
ggsave(sprintf("%s/figs/%s.d.png", dat1.dir, id), p, width = 8, height = 8)
