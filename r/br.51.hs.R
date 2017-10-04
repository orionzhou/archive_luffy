source("br.fun.R")

dirw = '/home/springer/zhoux379/scratch/briggs2/61.corncyc'

dirg = file.path(Sys.getenv("genome"), "Zmays_v4")
fg = file.path(dirg, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,16:18)]
gb = group_by(tg, par)
tg2 = summarise(gb, chr = chr[1], beg = min(beg), end = max(end), srd = srd[1], fam = names(sort(table(cat3), decreasing = T))[1])
tg2 = tg2[order(tg2$chr, tg2$beg),]

tg3 = tg2[tg2$chr == 1 & tg2$beg >= 25375937 & tg2$end <= 26296895,]
tg3$fam

tg3 = tg2[tg2$chr == 1 & tg2$beg >= 32272991 & tg2$end <= 34407017,]
tg3$fam

tx = merge(tm, tg3, by.x = 'gid', by.y = 'par')