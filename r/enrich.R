require(plyr)
require(tidyr)
require(dplyr)

fgo = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/15.tsv'
tgo = read.table(fgo, sep = "\t", as.is = T, header = F, quote = '')
colnames(tgo) = c("gid", "goid", "goname")

fg = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/gids.txt'
tg = read.table(fg)
gids.all = tg$V1


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

dirw = '/home/springer/zhoux379/data/misc2/peter_chg_go'

fi = file.path(dirw, '01.gid.txt')
ti = read.table(fi)
gids = unique(ti$V1)

fi = file.path(dirw, '02.gids.bg.txt')
ti = read.table(fi)
gids.bg = unique(ti$V1)


goid = "GO:0000723"

#study_size = length(gids)
#study_hit = sum(gids %in% tgo$gid[tgo$goid == goid])

study_size = length(gids.bg)
study_hit = sum(gids.bg %in% tgo$gid[tgo$goid == goid])

#pop_size = length(gids.bg)
#pop_hit = sum(gids.bg %in% tgo$gid[tgo$goid == goid])

pop_size = length(gids.all)
pop_hit = sum(tgo$goid == goid)


hitInSample = 8
hitInPop = 30 
failInPop = 39005 - hitInPop 
sampleSize = 4053


phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= F)

fisher.test(matrix(c(hitInSample, hitInPop-hitInSample, sampleSize-hitInSample, failInPop-sampleSize +hitInSample), 2, 2), alternative='greater')

pval.raw = phyper(study_hit-1, pop_hit, pop_size-pop_hit, study_size, lower.tail=F)

sprintf("study: %d/%d; population: %d/%d; fisher %g", study_hit, study_size, pop_hit, pop_size, pval.raw)