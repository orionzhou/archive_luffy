source("br.fun.R")

dirw = '/home/springer/zhoux379/data/misc2/grn23/47.coexp.test'

fi = file.path(dirw, "../37.rpkm.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

expr = t(as.matrix(ti[,-1]))
colnames(expr) = ti[,1]
gids = ti[,1]
ng = length(gids)

fc = file.path(dirw, "params.tsv")
tc = read.table(fc, sep = "\t", header = T, as.is = T)

for (i in 9:16) {

	nid = tc$nid[i]; coex_opt = tc$coex_opt[i]; n_edges = tc$n_edges[i]
	algorithm = tc$algorithm[i]; mcl_param = tc$mcl_param[i]
	
	net = sprintf("%s_%dk_%s", coex_opt, n_edges/1000, algorithm)
	if(algorithm == "mcl") {
		net = sprintf("%s_%g", net, mcl_param)
	}
	
	fd = sprintf("%s/01.edgeweight/%s.rda", dirw, coex_opt)
	if(i %in% c(1, 9, 17)) {
		x = load(fd)
		if(coex_opt == 'R') {
			ord = order(coexv)
		} else {
			ord = order(-coexv)
		}
	}
	idxs = ord[1:n_edges]
	dw = data.frame(g1 = gids[rowidx[idxs]], g2 = gids[colidx[idxs]], coex = coexv[idxs])
	#length(unique(c(dw$g1, dw$g2)))
	if(coex_opt == 'R') {
		dw$coex = exp((-(dw$coex-1)/25))
	}
	
	if(algorithm == 'mcl') {
		fo1 = sprintf("%s/05.modules/%s.1.tsv", dirw, net)
		write.table(dw, fo1, sep = "\t", row.names = F, col.names = F, quote = F)

		fo2 = sprintf("%s/05.modules/%s.2.txt", dirw, net)
		cmd = sprintf("mcl %s --abc -scheme 7 -I %g -o %s", fo1, mcl_param, fo2)
		system(cmd)
	
		fo3 = sprintf("%s/05.modules/%s.tsv", dirw, net)
		cmd = sprintf("mcl2tsv.py %s %s", fo2, fo3)
		system(cmd)
	} else {
		stopifnot(algorithm == 'clusterone')
		fo1 = sprintf("%s/05.modules/%s.1.txt", dirw, net)
		write.table(dw, fo1, sep = " ", row.names = F, col.names = F, quote = F)
	
		fo2 = sprintf("%s/05.modules/%s.2.csv", dirw, net)
		cmd = sprintf("java -jar $src/cluster_one-1.0.jar -f edge_list -F csv %s > %s", fo1, fo2)
		system(cmd)

		fo3 = sprintf("%s/05.modules/%s.tsv", dirw, net)
		cmd = sprintf("one2tsv.py %s %s", fo2, fo3)
		system(cmd)
	}
	
	fm = sprintf("%s/05.modules/%s.tsv", dirw, net)
	tm = read.table(fm, header = F, as.is = T, sep = "\t")
	cat(sprintf("%s: %d modules with %d genes\n", net, length(unique(tm$V1)), nrow(tm)))
}

for (i in 1:nrow(tc)) {
	nid = tc$nid[i]; coex_opt = tc$coex_opt[i]; n_edges = tc$n_edges[i]
	algorithm = tc$algorithm[i]; mcl_param = tc$mcl_param[i]
	
	net = sprintf("%s_%dk_%s", coex_opt, n_edges/1000, algorithm)
	if(algorithm == "mcl") {
		net = sprintf("%s_%g", net, mcl_param)
	}
	fm = sprintf("%s/05.modules/%s.tsv", dirw, net)
	tm = read.table(fm, header = F, as.is = T, sep = "\t")
	cat(sprintf("%s: %d modules with %d genes\n", net, length(unique(tm$V1)), nrow(tm)))
}