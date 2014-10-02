require(plyr)

dfam = file.path(Sys.getenv("data"), 'db', 'pfam')
ffam = file.path(dfam, 'genefam.tbl')
tfam = read.table(ffam, sep = "\t", header = T, as.is = T)

fg = file.path(Sys.getenv("genome"), 'HM101', 'augustus', '34.tbl')
tg = read.table(fg, sep = "\t", header = T, as.is = T)

tg = merge(tg, tfam, by.x = 'hid', by.y = 'dom', all.x = T)
tgs = tg[! tg$fam %in% unique(tfam$fam), ]
x = sort(table(tgs$hid), decreasing = T)
x[1:80]