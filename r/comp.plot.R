require(Gviz)
library(rtracklayer)
source("comp.plot.fun.R")

qname = "hm056"
tname = "hm101"

q = read_genome_stat(qname)
t = read_genome_stat(tname)
c = read_comp_stat(qname, tname)

tid = build_ideogram(t$gap, t$len)

f_mapp = "/home/youngn/zhoup/Data/genome/pan3/18_stat_k60/15_mapp.bw"

#### plot
fg = file.path(t$dir, "51.gtb")
tg = read.table(fg, header=T, sep="\t", quote="", as.is=T)[,c(1,3:6,15:17)]
tgs = tg[grepl('CRP', tg$cat3) | grepl('NBS', tg$cat3),]

for (i in 1:100) {
  chr = tgs$chr[i]
  beg = tgs$beg[i]
  end = tgs$end[i]
  fo = sprintf("%s/figs/genes/%s.png", c$dir, tgs$id[i])
  comp_plot(chr, beg, end, fo, q, t, c, tid, f_mapp)
}

########## comparison plot using a region in genome2
f_reg = '/home/youngn/zhoup/Data/genome/HM101/81_regions.tbl'
reg = read.table(f_reg, header=T, sep="\t", as.is=T)
i = 5
chr=reg$chr[i]; beg=reg$beg[i]; end=beg+30000


########## dotplot
xb = c(); xe = c(); yb = c(); ye = c()
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
tt = cbind(tbs, xb=xb, xe=xe, yb=yb, ye=ye)
p <- ggplot(tt) +
  geom_segment(mapping=aes(x=xb, xend=xe, y=yb, yend=ye), size=0.6) +
#  geom_segment(mapping=aes(x=xb, xend=xe, y=yb, yend=ye, color=factor(idc)), size=0.6) +
#  layer(data=dat2.coord, geom='rect', mapping=aes(xmin=hBeg.r, xmax=hEnd.r, ymin=-3000000, ymax=0, fill=hId), geom_params=list(size=0)) +
#  layer(data=dat2.coord, geom='text', mapping=aes(x=(hBeg.r+hEnd.r)/2, y=-4000000, label=hId), geom_params=list(hjust=0.5, vjust=1, angle=0, size=3)) +
  facet_grid(qId ~ hId, scales="free") +
#  scale_color_brewer(palette="Set1") +
  scale_x_continuous(name='HM101 (Mt4.0)') +
  scale_y_continuous(name=dat1.name) +
  theme_bw()
#  theme(legend.position="none")
#  theme(axis.text.x=element_blank(), axis.text.y=element_blank())
p
ggsave(sprintf("%s/figs/%s.d.png", dat1.dir, id), p, width=8, height=8)
