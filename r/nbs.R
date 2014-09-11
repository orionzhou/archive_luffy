dir = "/home/youngn/zhoup/Data/genome/HM101"
fp = file.path(dir, 'raw/34.pfam.tbl')
tp = read.table(fp, header = T, as.is = T, sep = "\t", quote = "")
tps = tp[tp$fam %in% c('PF00931'),] # , 'PF01582'

fj = file.path(dir, 'raw/26.longest.gtb')
tj = read.table(fj, header = T, as.is = T, sep = "\t", quote = "")[,c(1:6,18)]
fjs = file.path(dir, '42_nbs.gtb')
tjs = read.table(fjs, header = F, as.is = T, sep = "\t", quote = "")[,c(1:6,18)]
colnames(tjs) = c('id', 'par', 'chr', 'beg', 'end', 'srd', 'note')


ids_p = unique(tps$id)
ids_j = unique(tjs$id)

length(ids_p)
length(ids_j)
sum(ids_p %in% ids_j)

ids_pnj = ids_p[!ids_p %in% ids_j]
ids_jnp = ids_j[!ids_j %in% ids_p]
length(ids_pnj)
length(ids_jnp)

tp[tp$id %in% ids_jnp, ]
tj[tj$id %in% ids_jnp, ]

tp[tp$id %in% ids_pnj, ]
tj[tj$id %in% ids_pnj, ]

##### find all NBS-LRRs
ta = read.table(fj, header = T, as.is = T, sep = "\t", quote = "")
tn = ta[ta$id %in% c(ids_j, ids_pnj), ]

dir = "/home/youngn/zhoup/Data/misc2/nbs"
fo = file.path(dir, 'hm101.nbslrr.gtb')
write.table(tn, fo, sep = "\t", row.names = F, col.names = T, quote = F)

#####
dir = "/home/youngn/zhoup/Data/genome/HM101"
fj = file.path(dir, "raw/31.gtb")
tj = read.table(fj, header = T, as.is = T, sep = "\t", quote = "")[,c(1:6,18)]
fjs = file.path(dir, '42.nbs.gtb')
tjs = read.table(fjs, header = F, as.is = T, sep = "\t", quote = "")[,c(1:6,18)]
colnames(tjs) = c('id', 'par', 'chr', 'beg', 'end', 'srd', 'note')


fp = file.path(dir, '42.nbs/12.gtb')
tps = read.table(fp, header = T, as.is = T, sep = "\t", quote = "")[,c(1:6,18)]

ids_p = unique(tps$id)
ids_j = unique(tjs$id)
ids_pnj = ids_p[!ids_p %in% ids_j]
ids_jnp = ids_j[!ids_j %in% ids_p]

length(ids_p)
length(ids_j)
sum(ids_p %in% ids_j)
length(ids_pnj)
length(ids_jnp)


tj[tj$id %in% ids_jnp, ]
tj[tj$id %in% ids_pnj, ]

