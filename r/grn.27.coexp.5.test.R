source("br.fun.R")

dirw = '/home/springer/zhoux379/data/misc2/grn23/47.coexp.test'
#dirw = '/home/springer/zhoux379/scratch/briggs2/47.coexp.test'

fi = file.path(dirw, "../37.rpkm.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

expr = t(as.matrix(ti[,-1]))
colnames(expr) = ti[,1]
gids = ti[,1]
ng = length(gids)

## read GO
fg = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/15.tsv'
tg = read.table(fg, header = F, as.is = T, sep = "\t", quote = "")
colnames(tg) = c("gid", "goid", "goterm")
tg = tg[tg$gid %in% gids,]
tf = unique(tg[,c("gid","goterm")])
colnames(tf)[2] = "funcat"
tfs = ddply(tf, .(funcat), summarise, size = length(gid))
funcats = tfs$funcat[tfs$size >= 5]
t_go = tf[tf$funcat %in% funcats,]

## read GO - lev2
fg = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/15_lev2.tsv'
#tg = read.table(fg, header = F, as.is = T, sep = "\t")
#colnames(tg) = c("gid", "goid", "goterm")
#tg = tg[tg$gid %in% gids,]
#tf = unique(tg[,c("gid","goterm")])
#colnames(tf)[2] = "funcat"
#tfs = ddply(tf, .(funcat), summarise, size = length(gid))
#tf = tf[tf$funcat %in% tfs$funcat[tfs$size >= 10],]

## read CornCyc
fc = '/home/springer/zhoux379/data/genome/Zmays_v4/corncyc/01.tsv'
tc = read.table(fc, header = F, sep = "\t", quote = "")
gidlst = strsplit(tc$V7, split = ",")
lens = sapply(gidlst, length)
tc = data.frame(pl = rep(tc$V1, lens), 
	p1 = rep(tc$V2, lens), p2 = rep(tc$V3, lens), 
	p3 = rep(tc$V4, lens), p4 = rep(tc$V5, lens), 
	ps = rep(tc$V6, lens), gid = unlist(gidlst), stringsAsFactors = F)
table(tc$p1)

tp = tc[tc$p1 %in% c('Biosynthesis', 'Degradation/Utilization/Assimilation'),]
table(tp$p2)

tp2 = tp[tp$gid %in% gids, c('gid','p2')]
tf = unique(tp2)
colnames(tf) = c('gid','funcat')

tps = ddply(tc[tc$gid %in% gids,], .(pl), summarise, size = length(pl))
funcats = tps$pl[tps$size >= 5]
t_corncyc = tc[tc$pl %in% funcats, c('pl','gid')]
colnames(t_corncyc)[1] = 'funcat'

### Co-expression enrichment in GO categories / CornCyc Pathways (Z-test)
require(BSDA)

make_hist <- function(vec, xmin, xmax, xitv) {
	xb = seq(xmin, xmax, by = xitv)
	xm <- xb[-length(xb)] + 0.5 * diff(xb)
	td = data.frame(x = xm)
	
	x = findInterval(vec, xb)
	d = as.data.frame.table(table(x))
	d$x = as.numeric(levels(d$x))[d$x]
	d = d[d$x > 0 & d$x <= length(xm),]
	d$x = xm[d$x]
	d[!is.na(d$x),]
	td = merge(td, d, by = 'x', all.x = T)
	td[is.na(td[,2]), 2] = 0
	colnames(td)[2] = "freq"
	td
}

fun_sets = c("GO category", "CornCyc Pathway")

net = "C"
fd = sprintf("%s/01.edgeweight/%s.rda", dirw, net)
x = load(fd)

for (i in 1:2) {

	fun_set = fun_sets[i]
	label = strsplit(fun_set, split = " ")[[1]][1]
	tfun = tf
	if(label == "CornCyc") { tfun = tp }
	funcats = unique(tfun$funcat)

coexv_f = c()
coexv_r = c()
for (funcat in funcats) {
	gids_f = tfun$gid[tfun$funcat == funcat]
	idxs_f = which(gids %in% gids_f)
	idxs_r = sample(which(! gids %in% gids_f), length(idxs_f))
	coexm_f = coexm[idxs_f, idxs_f]
	coexm_r = coexm[idxs_r, idxs_r]
	coexv_f = c(coexv_f, coexm_f[lower.tri(coexm_f)])
	coexv_r = c(coexv_r, coexm_r[lower.tri(coexm_r)])
}

if (net == "C") {
	td1 = make_hist(coexv_f, -5, 5, 0.02)
	td2 = make_hist(coexv_r, -5, 5, 0.02)
	xlab = 'Co-expression Z-score'
}
if (net == "R") {
	coexv_f = exp((-(coexv_f-1)/20))
	coexv_r = exp((-(coexv_r-1)/20))
	td1 = make_hist(coexv_f, 0, 1, 0.002)
	td2 = make_hist(coexv_r, 0, 1, 0.002)
	xlab = "MutualRank Transformed Edge Weight"
}
if (net == "W") {
	td1 = make_hist(coexv_f, 0, 1, 0.002)
	td2 = make_hist(coexv_r, 0, 1, 0.002)
	xlab = "TOM (Topological Overlap Matrix) Score"
}

stopifnot(td1$x == td2$x)
td = cbind(td1, td2[,2])
colnames(td) = c('x', 'Funcat', 'Random')

p2 = ggplot(td) +
  geom_area(aes(x = x, y = Funcat, fill = 'Funcat'), alpha = 0.3) + 
  geom_area(aes(x = x, y = Random, fill = 'Random'), alpha = 0.3) + 
  #scale_x_continuous(name = xlab) +
  scale_x_continuous(name = xlab, limits = c(0.02,1)) +
  scale_y_continuous(name = 'Frequency') +
  scale_fill_manual(breaks = c("Funcat", "Random"), labels = c(fun_set, "Random"), values = c('red', 'blue')) +
  theme(legend.position = c(0.8,0.8), legend.key.size = unit(1, 'lines'), legend.title = element_blank(), legend.background = element_rect(fill = 'white', size=0)) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(axis.title.y = element_text(colour = 'pink', angle = 0)) +
  theme(axis.text.x = element_text(size = 9, color = "black")) +
  theme(axis.text.y = element_text(size = 9, color = "black"))

fp = sprintf("%s/11.goenrich/%s_%s.pdf", dirw, net, label)
ggsave(p2, filename = fp, width = 8, height = 8)

}

#### EGAD
require(EGAD)
annos = make_annotations(tf, gids, unique(tf$funcat))
	net = "N1"
	fed = sprintf("%s/02.edges/%s.txt", dirw, net)
	ted = read.table(fed, header = F, as.is = T, sep = " ")
	
	nw = build_binary_network(as.matrix(ted[,1:2]), gids)
	nd = node_degree(nw)
	gba_auc_nv = run_GBA(nw, annos, nFold = 3)
	aucA <- gba_auc_nv[,1]
	aucB <- gba_auc_nv[,3]
	plot_distribution(aucA)
	plot_density_compare(aucA, aucB)
	plot_value_compare(aucA, aucB)
	
	fp = sprintf("%s/11.goenrich/%s.pdf", dirw, net)
	pdf(fp, width = 6, height = 6)
	plot_distribution(res, xlab="Neighbor voting AUROC",
		ylab="Number of functional terms", 
		b=30, xlim=c(0.4,1), ylim=c(0, 440), col="gray64", density=F, avg=F, bar = T)
	dev.off()

#### enrichment in functional categories
total_gic = length(unique(tf$gid))
total_gnic = length(gids) - total_gic

nets = c(sprintf("C%d", 1:5), sprintf("N%d", 1:5))

for (net in nets) {
	fm = sprintf("%s/05.modules/%s.tsv", dirw, net)
	tm = read.table(fm, header = T, as.is = T, sep = "\t")

	tms = ddply(tm, .(grp), summarise, size = length(gid))
	sizes = tms$size
	cat(sprintf("%s: %g modules, %g genes, mean size %g +- %.01f, median size %g +- %.01f\n",
		net, nrow(tms), length(unique(tm$gid)), round(mean(sizes)), sd(sizes),
		round(median(sizes)), mad(sizes)))

	total_gim = length(unique(tm$gid))
	total_gnim = length(gids) - total_gim

tz = merge(tf, tm, by = 'gid', all.x = T)
tzs = ddply(tz, .(funcat), summarise, psize = length(gid), gim = sum(!is.na(grp)), nmodule = length(unique(grp[!is.na(grp)])), gnim = sum(is.na(grp)), pct_im = gim/total_gim*100, pct_nim = gnim/total_gnim*100)
tzs = tzs[tzs$psize > 10,]

pvals = apply(tzs, 1, myfunc <- function(x) dhyper(as.numeric(x[3]), total_gim - as.numeric(x[3]), total_gnim - as.numeric(x[5]), as.numeric(x[2])))
tzs = cbind(tzs, pval = p.adjust(pvals, 'BH'))
tzs = tzs[order(tzs$p2), c('p2','gim','psize','pct_im','pct_nim','pval')]
tzs = cbind(x=1:nrow(tzs), tzs)

tsig = tzs[tzs$pval < 0.05,]
tsig = cbind(tsig, xmin = tsig$x-0.25, xmax = tsig$x+0.25, y = pmax(tsig$pct_im, tsig$pct_nim)+0.2, lab = "*  ")
tsig$lab[tsig$pval < 0.005] = "** "
tsig$lab[tsig$pval < 0.0005] = "***"

to = gather(tzs[,c('funcat','pct_im','pct_nim')], type, percent, -funcat)
p1 = ggplot(to, aes(x = funcat, y = percent, fill = type)) +
  geom_bar(stat = 'identity', position = 'dodge', width = 0.7) +
  geom_signif(y_position=tsig$y, xmin=tsig$xmin, xmax=tsig$xmax, annotation=tsig$lab,
  	tip_length = 0.02, vjust = 2, hjust = -0.7) +
  #scale_x_continuous(name = '') +
  scale_y_continuous(expand = c(0.08, 0)) +
  scale_fill_manual(values = c("firebrick", "grey"), labels = c("in modules", "not in modules"), name = net) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = c(0.55, 0.04), legend.direction = "vertical", legend.justification = c(0,0), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 9), legend.box.background = element_rect(color='black')) +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = sprintf("%s/20.%s.pdf", dirw, net)
ggsave(p1, filename = fp, width = 6, height = 6)

}

## GO enrichment
fc = file.path(dirw, "params.tsv")
tc = read.table(fc, sep = "\t", header = T, as.is = T)

test_funcat_enrich <- function(rw, ng) {
	sizeo =  as.integer(rw[3]); sizem = as.integer(rw[4]); sizef = as.integer(rw[5])
		n1 = sizeo
		n2 = sizef
		n3 = ng - n2
		n4 = sizem
	dhyper(n1, n2, n3, n4)
}

for (i in 12:nrow(tc)) {
	nid = tc$nid[i]; coex_opt = tc$coex_opt[i]; n_edges = tc$n_edges[i]
	algorithm = tc$algorithm[i]; mcl_param = tc$mcl_param[i]
	
	net = sprintf("%s_%dk_%s", coex_opt, n_edges/1000, algorithm)
	if(algorithm == "mcl") {
		net = sprintf("%s_%g", net, mcl_param)
	}
	fm = sprintf("%s/05.modules/%s.tsv", dirw, net)
	tm = read.table(fm, header = T, as.is = T, sep = "\t")
	colnames(tm) = c("grp", "gid")
	tms = ddply(tm, .(grp), summarise, sizem = length(gid))
	sizes = tms$size
	
	tt = t_go
	tts = ddply(tt, .(funcat), summarise, sizef = length(gid))

	mods = tms$grp
	cats = tts$funcat
	nmod = length(mods)
	nfun = length(cats)
	to =  data.frame(mod = rep(mods, each = nfun), funcat = rep(cats, nmod), stringsAsFactors = F)
	
	to1 = merge(tm, tt, by = 'gid')
	to2 = ddply(to1, .(grp, funcat), summarise, sizeo = length(gid))
	colnames(to2)[1] = 'mod'
	to3 = merge(to, to2, by = c('mod','funcat'), all.x = T)
	to3$sizeo[is.na(to3$sizeo)] = 0
	to4 = merge(to3, tms, by.x = 'mod', by.y = 'grp')
	to5 = merge(to4, tts, by = 'funcat')
	stopifnot(nrow(to5) == nrow(to))
	to6 = to5[,c(2,1,3:5)]
	
	pvals = apply(to6, 1, test_funcat_enrich, ng)
	tv = cbind(to6, pval.adj = p.adjust(pvals, method = "BH"))
	tv2 = tv[tv$pval.adj < 0.05,]
	nm1 = length(unique(tv2$mod))
	nf1 = nrow(unique(tv2[,c(1:2)]))
	
	
	tt = t_corncyc
	tts = ddply(tt, .(funcat), summarise, sizef = length(gid))

	mods = tms$grp
	cats = tts$funcat
	nmod = length(mods)
	nfun = length(cats)
	to =  data.frame(mod = rep(mods, each = nfun), funcat = rep(cats, nmod), stringsAsFactors = F)
	
	to1 = merge(tm, tt, by = 'gid')
	to2 = ddply(to1, .(grp, funcat), summarise, sizeo = length(gid))
	colnames(to2)[1] = 'mod'
	to3 = merge(to, to2, by = c('mod','funcat'), all.x = T)
	to3$sizeo[is.na(to3$sizeo)] = 0
	to4 = merge(to3, tms, by.x = 'mod', by.y = 'grp')
	to5 = merge(to4, tts, by = 'funcat')
	stopifnot(nrow(to5) == nrow(to))
	to6 = to5[,c(2,1,3:5)]
	
	pvals = apply(to6, 1, test_funcat_enrich, ng)
	tv = cbind(to6, pval.adj = p.adjust(pvals, method = "BH"))
	tv2 = tv[tv$pval.adj < 0.05,]
	nm2 = length(unique(tv2$mod))
	nf2 = nrow(unique(tv2[,c(1:2)]))
	
	cat(sprintf("%s %g %g %g+%.01f %g+%.01f %d %d %d %d\n",
		net, nrow(tms), length(unique(tm$gid)), 
		round(mean(sizes)), sd(sizes),
		round(median(sizes)), mad(sizes),
		nm1, nf1, nm2, nf2
		))
}