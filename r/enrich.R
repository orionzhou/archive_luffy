require(plyr)
require(tidyr)
require(dplyr)

fg = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/15.tsv'
tg = read.table(fg, sep = "\t", as.is = T, header = F, quote = '')
colnames(tg) = c("gid", "goid")

fd = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/16.go.tsv'
td = read.table(fd, sep = "\t", header = T, as.is = T, quote = '')

fi = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/gids.txt'
ti = read.table(fi, as.is = T)
gids.all = ti$V1


go_enrich <- function(gids, gids.all. = gids.all, tg. = tg, td. = td) {
	tgs = tg %>% 
		group_by(goid) %>% 
		summarise(hitInPop = n())

	sampleSize = length(gids)

	tz = tg[tg$gid %in% gids,] %>%
		group_by(goid) %>% 
		summarise(hitInSample = n())

	tw = merge(tz, tgs, by = 'goid')
	tw = tw %>% 
		mutate(sampleSize = length(gids), popSize = length(gids.all), 
			pval.raw = phyper(hitInSample-1, hitInPop, popSize-hitInPop, sampleSize, lower.tail = F),
			pval.adj = p.adjust(pval.raw, method = "BH")) %>%
		filter(pval.raw < 0.05) %>%
		arrange(pval.adj) %>%
		transmute(goid = goid, ratioInSample = sprintf("%d/%d", hitInSample, sampleSize),
			ratioInPop = sprintf("%d/%d", hitInPop, popSize),
			pval.raw = pval.raw, pval.adj = pval.adj)

	to = merge(tw, td, by = 'goid')
	to %>% arrange(pval.adj, pval.raw)
}

#fisher.test(matrix(c(hitInSample, hitInPop-hitInSample, sampleSize-hitInSample, failInPop-sampleSize +hitInSample), 2, 2), alternative='two.sided')



dirw = '/home/springer/zhoux379/data/misc2/peter_chg_go'

fi = file.path(dirw, '01.gid.txt')
ti = read.table(fi, as.is = T)
gids = unique(ti$V1)

fi = file.path(dirw, '02.gids.bg.txt')
ti = read.table(fi, as.is = T)
gids = unique(ti$V1)

to = go_enrich(gids)


