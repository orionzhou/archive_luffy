require(plyr)
require(ggplot2)
require(dplyr)

dirw = '/home/springer/zhoux379/data/misc1/jackie.acr'

fb = file.path(dirw, "01.bed")
tb = read.table(fb, sep = "\t", header = T, as.is = T)
tb$start = tb$start + 1
tb2 = cbind.data.frame(tb, lid=sprintf("%d-%d-%d",tb$chr,tb$start,tb$stop), stringsAsFactors = F)
ids_all = unique(tb2$lid)

# cd $misc1/jackie.acr
# seqret.pl -d $genome/Zmays_v4/11_genome.fas -b 01.bed -o 02.fas
# blat $genome/PH207/21.blat/db.2bit 02.fas -ooc=$genome/PH207/21.blat/db.2bit.tile11.ooc 04.psl
# psl2tsv.pl -i 04.psl -o 05.tsv

fi = file.path(dirw, "05.tsv")
tir = read.table(fi, sep = "\t", header = T, as.is = T)
ti = tir[tir$alnLen/tir$qSize >= 0.8 & tir$ident >= 0.8,]

ids_mapped = unique(ti$qId)

grp = dplyr::group_by(ti, qId)
ti2 = dplyr::summarise(grp, nbest = sum(score == max(score)), score = max(score))
ti3 = merge(ti, ti2, by = c("qId", "score"))
ids_mapped_u = unique(ti2$qId[ti2$nbest == 1])
ids_mapped_m = unique(ti2$qId[ti2$nbest > 1])


# unique hits
tk = ti3[ti3$qId %in% ids_mapped_u,]
stopifnot(nrow(tk) == length(ids_mapped_u))
ids_mapped_u1 = tk$qId[tk$alnLen == tk$qSize & tk$misMatch == 0 & tk$qNumIns+tk$tNumIns==0]
ids_mapped_u2 = tk$qId[tk$alnLen < tk$qSize & tk$misMatch == 0 & tk$qNumIns+tk$tNumIns==0]
ids_mapped_u3 = tk$qId[tk$misMatch > 0 | tk$qNumIns+tk$tNumIns > 0]
stopifnot(length(ids_mapped_u1) + length(ids_mapped_u2) + length(ids_mapped_u3) == length(ids_mapped_u))

tm = tk[tk$qId %in% ids_mapped_u3,]
ids_mapped_u31 = tm$qId[tm$misMatch > 0 & tm$qNumIns+tm$tNumIns == 0]
ids_mapped_u32 = tm$qId[tm$misMatch == 0 & tm$qNumIns+tm$tNumIns > 0]
ids_mapped_u33 = tm$qId[tm$misMatch > 0 & tm$qNumIns+tm$tNumIns > 0]
stopifnot(length(ids_mapped_u31) + length(ids_mapped_u32) + length(ids_mapped_u33) == length(ids_mapped_u3))


outputs = c(
'',
sprintf("%5d total sequences:", length(ids_all)),
sprintf("%5d do not map to PH207 (0.8 identity, 0.8 coverage)", length(ids_all)-length(ids_mapped)),
sprintf("%5d mapped >=1 times to PH207:", length(ids_mapped)),
sprintf("   %5d mapped uniquely:", length(ids_mapped_u)),
sprintf("   %5d mapped multiple times:", length(ids_mapped_m)),
''
)
cat(paste(outputs, collapse = "\n"))


### read synteny info
fg = '/home/springer/zhoux379/data/genome/Zmays_v4/51.gtb'
tg = read.table(fg, sep = "\t", header = T, as.is = T)


fm = '/home/springer/zhoux379/data/genome/PH207/synteny/mapping.tsv'
tm = read.table(fm, sep = "\t", header = F, as.is = T)
colnames(tm) = c("bchr","bbeg",'bend','bgid','pchr','pbeg','pend','pgid','type1','type2')
tm2 = tm[,c('bgid','pgid','pchr','pbeg','pend','type1','type2')]
idxs = which(tm2$pchr %in% sprintf("%d", 1:10))
tm2$pchr[idxs] = sprintf("chr%02d", as.integer(tm2$pchr[idxs]))
tm3 = tm2[tm2$type2 %in% c('One-to-One'),]

tj = ti3
tj2 = merge(tj, tb2, by.x = 'qId', by.y = 'lid')
stopifnot(nrow(tj2) == nrow(tj))

length(unique(tj2$gene))
sum(unique(tj2$gene) %in% tm3$bgid)
tj4 = merge(tj2, tm3, by.x='gene', by.y='bgid', all.x = T)
stopifnot(nrow(tj4) == nrow(tj2))
tj4 = within(tj4, {
    pdist = ifelse(tId == pchr, 
        ifelse(tEnd < pbeg, pbeg - tEnd,
            ifelse(tBeg > pend, tBeg - pend, 0)), -1)
})
tj4$pdist[is.na(tj4$pdist)] = -1
tj4 = within(tj4, {
    syntag = ifelse(pdist >= 0 & pdist <= 1000000, 'inSynteny', 'notInSynteny')
})


tj5 = tj4[tj4$qId %in% ids_mapped_m,]
ids_mapped_m_resolved = unique(tj5$qId[tj5$syntag])
ids_mapped_m_unresolved = ids_mapped_m[! ids_mapped_m %in% ids_mapped_m_resolved]
tj6 = tj5[!(tj5$qId %in% ids_mapped_m_resolved & tj5$syntag == 'notInSynteny'), ]
#pick first occurrence of each qId
tj7 = tj6[match(unique(tj6$qId), tj6$qId), ]


tj8 = rbind(tj4[tj4$qId %in% ids_mapped_u,], tj7)
tj8 = within(tj8, {
	vartag = ifelse(alnLen == qSize & misMatch == 0 & qNumIns+tNumIns == 0, 'Identical', 
	ifelse(alnLen < qSize & misMatch == 0 & qNumIns+tNumIns == 0, 'assemblyGap',
	ifelse(misMatch > 0 & qNumIns+tNumIns == 0, 'onlySNP', 'Indel')))
	maptag = ifelse(nbest == 1, 'Unique', "multiMapping")
})

to = tj8
outputs = c(
'',
sprintf("%5d total sequences:", length(ids_all)),
sprintf("%5d do not map to PH207 (0.8 identity, 0.8 coverage)", length(ids_all)-length(ids_mapped)),
sprintf("%5d mapped >=1 times to PH207:", nrow(to)),
sprintf("   %5d mapped uniquely", sum(to$maptag=='Unique')),
sprintf("   %5d mapped multiple times", sum(to$maptag=='multiMapping')),
'',
sprintf("%5d mapped >=1 times to PH207:", nrow(to)),
sprintf("   %5d have the mapped location near PH207 ortholog", sum(to$syntag=='inSynteny')),
sprintf("   %5d don't have a PH207 ortholog or are mapped to a distant location from the PH207 ortholog", sum(to$syntag=='notInSynteny')),
'',
sprintf("%5d mapped >=1 times to PH207:", nrow(to)),
sprintf("   %5d full conservation (100%% coverage, no SNP, no InDel)", sum(to$vartag=='Identical')),
sprintf("   %5d likely near assembly gaps (100%% coverege, no SNP, no InDel)", sum(to$vartag=='assemblyGap')),
sprintf("   %5d only SNP(s), no Indel", sum(to$vartag=='onlySNP')),
sprintf("   %5d >=1 InDel(s) (and/or SNPs)", sum(to$vartag=='Indel')),
''
)
cat(paste(outputs, collapse = "\n"))


outputs = c(
'',
sprintf("Of the %5d ACRs with at least one InDel(s) (and/or SNPs):", sum(to$vartag=='Indel')),
sprintf("   %5d have >= 100 inserted B73 bases (deleted PH207 bases)", sum(to$vartag=='Indel' & to$qBaseIns >= 100)),
sprintf("       %5d exactly one B73 insertion event (PH207 deletion)", sum(to$vartag=='Indel' & to$qBaseIns >= 100 & to$qNumIns == 1)),
sprintf("   %5d have >= 100 deleted B73 bases (inserted PH207 bases)", sum(to$vartag=='Indel' & to$tBaseIns >= 100)),
sprintf("       %5d exactly one B73 deletion event (PH207 insertion)", sum(to$vartag=='Indel' & to$tBaseIns >= 100 & to$tNumIns == 1)),

''
)
cat(paste(outputs, collapse = "\n"))

tp = to[,c('chr','start','stop','summit','gene','distance','category', 
	'qId','qBeg','qEnd','qSrd','qSize',
	'tId','tBeg','tEnd','tSrd','tSize',
	'alnLen','match','misMatch','baseN','qNumIns','tNumIns','qBaseIns','tBaseIns','ident',
	'nbest','pgid','pdist','syntag','vartag','maptag')]
colnames(tp)[colnames(tp) == 'pgid'] = 'gene.PH207'
colnames(tp)[colnames(tp) == 'pdist'] = 'distance.PH207'
fo = file.path(dirw, "09.tsv")
write.table(tp, fo, sep = "\t", row.names = F, col.names = T, quote = F)


