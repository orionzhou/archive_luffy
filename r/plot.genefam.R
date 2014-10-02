require(ggplot2)
require(RColorBrewer)
source("loc.R")

ff = '/home/youngn/zhoup/Data/misc2/genefam/nbs.info'
tf = read.table(ff, sep="\t", header=T, as.is=T)[,1:2]
labs = unique(tf$id)
length(labs)

ocols = brewer.pal(6, "Dark2")
oshas = 0:5 #15:19

cols = rep(ocols, 6)
shas = rep(oshas, each = 6)

##### NBS-LRR
## HM101
org = "HM101"
dir = file.path('/home/youngn/zhoup/Data/genome', org, '42.nbs')
fi = file.path(dir, "12.gtb")
t1 = read.table(fi, header = T, as.is = T, sep = "\t", quote = "")[,c(1,3:5,17)]
t1 = cbind(t1, rpos = (t1$beg + t1$end) / 2)

##non-HM101
org = "HM018"
dir = file.path('/home/youngn/zhoup/Data/genome', org, '42.nbs')
fi = file.path(dir, "16.liftover.tbl")
t1 = read.table(fi, header = T, as.is = T, sep = "\t", quote = "")[,c(1,6,9:10)]
colnames(t1) = c('id', 'cat3', 'chr', 'rpos')

## plot
t1$cat3 = factor(t1$cat3, levels = labs)

chrs = sort(unique(t1$chr))
tloc = data.frame(chr = chrs, chrnum = 1:length(chrs))

t2 = merge(t1, tloc, by = 'chr')
locs = t2$chrnum * 1000000000 + t2$rpos
names(locs) = t2$id
cl = locCluster(locs, 50000)
tlo = merge(t2, cl, by="id")

p <- ggplot(data = tlo) +
  geom_point(mapping = aes(x = rpos/1000000, y = cluster_y,
    colour = cat3, shape = cat3), size = 0.9) +
  scale_x_continuous(name='Chr Position (Mbp)', expand=c(0.01, 0)) +
  scale_y_continuous(name='', expand=c(0.04, 0)) +
  facet_grid(chr ~ .) + 
  scale_colour_manual(name = "cat3", breaks = labs, labels = labs, values = cols) +
  scale_shape_manual(name = "cat3", breaks = labs, labels = labs, values = shas) +
  theme_bw() +
  theme(legend.position='right', legend.title = element_blank()) +
  guides(color = guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(size=8, angle=0)) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank())
do = '/home/youngn/zhoup/Data/misc2/genefam/fig.nbs'
ggsave(p, filename=sprintf("%s/%s.pdf", do, org), width=15, heigh=9)

