source("comp.plot.fun.R")

########## load genome1 datasets
dat1.name = "hm056"
#dat1.name = "hm340"

dat1.dir = file.path("/home/youngn/zhoup/Data/misc3", dat1.name)

# sequence lengths
f_len = file.path(dat1.dir, "11_seqlen.tbl")
dat1.seqlen = read.table(f_len, header=TRUE, sep="\t", as.is=T)

# gap locations
f_gap = file.path(dat1.dir, "12_gaploc.tbl")
dat1.gap = read.table(f_gap, header=TRUE, sep="\t", as.is=T)

# gene annotation
dat1.gene = NULL
dat1.te = NULL
dat1.nbs = NULL

#f_crp = "/home/youngn/zhoup/Data/misc4/spada.crp.HM340/31_model_evaluation/61_final.gtb"
#dat1.crp = read.table(f_crp, header=TRUE, sep="\t", as.is=T, quote=NULL)[,c(3:6,18)]
dat1.crp = NULL

# mappability
dat1.mapp = NULL

# wrapping
dat1 = list( name=dat1.name, dir=dat1.dir, seqlen=dat1.seqlen, gap=dat1.gap,
	gene=dat1.gene, te=dat1.te, nbs=dat1.nbs, crp=dat1.crp, mapp=dat1.mapp )

########## load genome2 datasets
dat2.name = "hm101"

dat2.dir = "/home/youngn/zhoup/Data/genome/Mtruncatula_4.0"

# sequence lengths
f_len = file.path(dat2.dir, "15_seqlen.tbl")
dat2.seqlen = read.table(f_len, header=TRUE, sep="\t", as.is=T)

# gap locations
f_gap = file.path(dat2.dir, "16_gap_loc.tbl")
dat2.gap = read.table(f_gap, header=TRUE, sep="\t", as.is=T)

# gene annotation
f_gen = file.path(dat2.dir, "21_gene.gtb")
dat2.gene = read.table(f_gen, header=TRUE, sep="\t", as.is=T, quote=NULL)[,c(3:6,18)]

f_te = file.path(dat2.dir, "41_te.gtb")
dat2.te = read.table(f_te, header=TRUE, sep="\t", as.is=T, quote=NULL)[,c(3:6,18)]

f_nbs = file.path(dat2.dir, "42_nbs.gtb")
dat2.nbs = read.table(f_nbs, header=TRUE, sep="\t", as.is=T, quote=NULL)[,c(3:6,18)]

f_crp = file.path(dat2.dir, "43_crp.gtb")
dat2.crp = read.table(f_crp, header=TRUE, sep="\t", as.is=T, quote=NULL)[,c(3:6,18)]

# mappability
dat2.mapp = NULL #import("/home/youngn/zhoup/Data/db/gem/Mtruncatula_4.0_90mer.bw")

# wrapping
dat2 = list( name=dat2.name, dir=dat2.dir, seqlen=dat2.seqlen, gap=dat2.gap,
	gene=dat2.gene, te=dat2.te, nbs=dat2.nbs, crp=dat2.crp, mapp=dat2.mapp )

########## load (blast) comparison dataset
tc1 = read.table(file.path(dat1.dir, "23_blat/17_chain.gal"), header=TRUE, sep="\t", as.is=T)
tb1 = read.table(file.path(dat1.dir, "23_blat/18_block.gal"), header=TRUE, sep="\t", as.is=T)

tc=tc1
tb=tb1

########## comparison plot using a region in genome2
regions = read.table(file.path(dat2.dir, "81_regions.tbl"), header=T, sep="\t", as.is=T)
i = 8
chr=regions$chr[i]; beg=regions$beg[i]; end=regions$end[i]

# save path for PNG file
fn = sprintf("%s/figs/%s_%d.png", dat1.dir, chr, beg)

# filtering comparison
tas = aln[aln$hId==chr & ( (beg<=aln$hBeg & aln$hBeg<=end) | (beg<=aln$hEnd & aln$hEnd<=end) ), ]
source("comp.plot.fun.R")

# data preprocessing
dat = data_preprocess(tas, dat1, dat2)

# plot
source("comp.plot.fun.R")
subtitle = sprintf("%s:%g-%gK", chr, beg/1000, end/1000)
comp.plot(fn, dat, width=2000, height=1000, subtitle=subtitle)

########## batch plot
regions = read.table(file.path(dat2.dir, "81_regions.tbl"), header=T, sep="\t", as.is=T)
for (i in 1:nrow(regions)) {
	chr=regions$chr[i]; beg=regions$beg[i]; end=regions$end[i]
	fn = sprintf("%s/figs/%s_%d.png", dat1.dir, chr, beg)
	tas = aln[aln$hId==chr & ( (beg<=aln$hBeg & aln$hBeg<=end) | (beg<=aln$hEnd & aln$hEnd<=end) ), ]
	dat = data_preprocess(tas, dat1, dat2)
	subtitle = sprintf("%s:%g-%gK", chr, beg/1000, end/1000)
	comp.plot(fn, dat, width=2000, height=1000, subtitle=subtitle)
}

########## comparison plot using a region in genome1
id = "scaffold_314"
tcs = tc[tc$qId==id,]
tbs = tb[tb$qId==id,]

	dat = data_preprocess(tcs, tbs, dat1, dat2)
	subtitle = sprintf("%s", id)
	fn = sprintf("%s/figs/%s.png", dat1.dir, id)
	comp.plot(fn, dat, width=2000, height=1000, subtitle=subtitle)

fcn = file.path(dat1.dir, "23_blat/62_chain.gal")
tcn = read.table(fcn, header=TRUE, sep="\t", as.is=T, quote=NULL)
fbn = file.path(dat1.dir, "23_blat/63_block.gal")
tbn = read.table(fbn, header=TRUE, sep="\t", as.is=T, quote=NULL)
tcs = tcn[tcn$qId == id,]
tbs = tbn[tbn$qId == id,]

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
