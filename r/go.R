require(dplyr)

fgo = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/15.tsv'
tgo = read.table(fgo, sep = "\t", as.is = T, header = F, quote = '')
colnames(tgo) = c("gid", "goid", "goname")

go_enrich <- function(gids) {
	grp = dplyr::group_by(tgo, goid)
	tgos = dplyr::summarise(grp, size_go = n(), goname=goname[1])

	size_gene = sum(gids %in% tgo$gid)

	tz1 = tgo[tgo$gid %in% gids,]
	grp = dplyr::group_by(tz1, goid)
	tz2 = dplyr::summarise(grp, sizeo = n())
	tz3 = data.frame(goid = tgos$goid[!tgos$goid %in% tz2$goid], sizeo = 0, stringsAsFactors = F)
	tz4 = rbind(data.frame(tz2), tz3)

	tz4 = merge(tz4, tgos, by = 'goid')
	tz4 = cbind(tz4, size_gene = size_gene, total_gene = length(unique(tgo$gid)))

	tz5 = within(tz4, {
		pval.raw = phyper(sizeo-1, size_go, total_gene-size_go, size_gene, lower.tail=F)
	})
	tz6 = cbind(tz5, pval.adj = p.adjust(tz5$pval.raw, method = "BH"))
	tz = tz6[tz6$pval.adj <= 0.05,c('goid','sizeo','size_go','pval.raw','pval.adj','goname')]
	tz = tz[order(tz$pval.adj, decreasing = F),]
	tz
}

