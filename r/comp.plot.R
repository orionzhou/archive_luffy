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
f_mapp = "/home/youngn/zhoup/Data/db/gem/Mtruncatula_4.0.bw"
dat2.mapp = import(f_mapp)

# wrapping
dat2 = list( name=dat2.name, dir=dat2.dir, seqlen=dat2.seqlen, gap=dat2.gap,
	gene=dat2.gene, te=dat2.te, nbs=dat2.nbs, crp=dat2.crp, mapp=dat2.mapp )

########## load (blast) comparison dataset
f_aln = file.path(dat1.dir, '21_blastn/05_tiled.tbl')
aln = read.table(f_aln, header=TRUE, sep="\t", as.is=T)


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
id = "scaffold_119"

tas = aln[aln$qId==id & aln$qLen >= 1000,]

