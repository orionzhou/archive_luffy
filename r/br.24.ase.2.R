require(tidyverse)
require(RColorBrewer)

dirw = file.path(Sys.getenv("misc2"), "briggs", "46.ase")

fm = file.path(dirw, '../00.1.read.tsv')
tm = read.table(fm, header = T, sep = "\t", as.is = T)[,1:5]


## Indel mapping bias
sids = c("BR006", "BR054", "BR151")
opts = c("snp+indel", "snp")

for (sid in sids) {
for (opt in opts) {


fi = sprintf("%s/%s.5.%s.bed", diri, sid, opt)
ti = read.table(fi, header = F, sep = "\t", as.is = T)
colnames(ti) = c("chr", "beg", "end", "rid", "nref", "nalt", "nunk", "nerr")
ti = ti[ti$nref+ti$nalt>0, ]

n0 = nrow(ti)

n1 = sum(ti$nref + ti$nalt == 1)
n11 = sum(ti$nref == 1 & ti$nalt == 0)
n12 = sum(ti$nref == 0 & ti$nalt == 1)

n2 = sum(ti$nref + ti$nalt >= 2)
n21 = sum(ti$nref > 1 & ti$nalt == 0)
n22 = sum(ti$nref == 0 & ti$nalt > 1)
n23 = sum(ti$nref > 0 & ti$nalt > 0)


lab = tl$label[tl$SampleID==sid]
outputs = c(
'',
sprintf("%s [%s] against %s", sid, lab, opt),
sprintf("%8d total reads with >=1 variants:", n0),
sprintf("%8d: only 1 variants", n1),
sprintf("   %8d [%5.02f%%] B73 allele", n11, n11/n1*100),
sprintf("   %8d [%5.02f%%] Mo17 allele", n12, n12/n1*100),
sprintf("%8d: >=2 variants", n2),
sprintf("   %8d [%5.02f%%] with all B73 allele", n21, n21/n2*100),
sprintf("   %8d [%5.02f%%] with all Mo17 allele", n22, n22/n2*100),
sprintf("   %8d [%5.02f%%] with both B and M alleles", n23, n23/n2*100),
''
)
cat(paste(outputs, collapse = "\n"))


}
}

## calc gene stats (should be in rnaseq.ase.py)
diri = "/home/springer/zhoux379/scratch/briggs/28.ase"
for (i in 166:219) {
	sid = tm$SampleID[i]
	f1 = sprintf("%s/%s.tsv", diri, sid)
	f2 = sprintf("%s/%s.bed", diri, sid)
	if(! file.exists(f1) | !file.exists(f2)) {
		cat(sid, 'skipped\n')
		next
	}
	fo = sprintf("%s/01.all/%s.tsv", dirw, sid)
	cmd = sprintf("bed.ase.sum.py %s %s %s", f1, f2, fo)
	system(cmd)
	cat(sid, "\n")
}

## genotyping
tm2 = cbind(tm, propb = NA)
for (i in 1:nrow(tm)) {
	sid = tm$SampleID[i]
	fi = sprintf("%s/01.all/%s.tsv", dirw, sid)
	ti = read_tsv(fi, col_names = T)
	ti = ti %>%
		filter(nref + nalt >= 30) %>%
		mutate(propb = nref/(nref+nalt))
	tm2$propb[i] = median(ti$propb)
	cat(sid, "\n")
}
tm3 = tm2 %>%
	mutate(gt = ifelse(propb <= 0.05, "Mo17", 
		ifelse(propb >= 0.95, "B73",
		ifelse(propb >= 0.45 & propb <= 0.55, "B73xMo17",
		ifelse(propb >= 0.30 & propb <= 0.40, "Mo17xB73",
		ifelse(propb >= 0.65 & propb <= 0.75, "B73xMo17", 'unknown'))))))
tm3$gt[is.na(tm3$gt)] = 'unknown'
tm3 = mutate(tm3, gt_ok = ifelse(gt == Genotype, '1', '3'))
table(tm3$gt, useNA='ifany')

fo = file.path(dirw, "13.gt.tsv")
write.table(tm3[,c('SampleID','propb','gt','gt_ok')], fo, sep = "\t", row.names = F, col.names = T, quote = F)


## concatenate results
tl1 = tl[tl$Genotype %in% c("B73xMo17", "Mo17xB73"),]

to = data.frame()
for (i in 1:nrow(tl1)) {
	sid = tl1$SampleID[i]
	geno = tl1$Genotype[i]
	tiss = tl1$Tissue[i]
	repl = tl1$Treatment[i]
	
	fi = sprintf("%s/01.all/%s.tsv", diro, sid)
	if(! file.exists(fi)) {
		cat(sid, 'skipped\n')
		next
	}
	ti = read.table(fi, sep = "\t", header = T, as.is = T)
	#ti = ti[ti$nref + ti$nalt >= 30, ]
	
	geno = ifelse(geno == "Mo17xB73", "MxB", "BxM")
	lab = sprintf("%s %s [%d]", tiss, geno, nrow(ti))
	
	ti = cbind(sid = sid, ti)
	to = rbind(to, ti)
	cat(sid, "\n")
}
fo = file.path(dirw, "10.ase.gene.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)


# plot ASE summary
fi = file.path(dirw, "10.ase.gene.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

p1 = ggplot(ti[ti$nref + ti$nalt >= 30, ]) +
  geom_histogram(aes(x = nref/(nref+nalt)), bins = 50) +
  scale_x_continuous(name = 'Proportion B73-specific Reads', limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous(name = "# Genes") +
  #scale_color_manual(name = "", values = cols) +
  facet_wrap( ~ sid, nrow = 8) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(legend.position = c(0.4, 0.8), legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = file.path(dirw, "11.ase.gene.pdf")
ggsave(p1, filename = fp, width = 13, height = 13)

# tissue level summary
fi = file.path(dirw, "10.ase.gene.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

tb = merge(ti, tl[,c("SampleID", "Tissue", "Genotype")], by.x = 'sid', by.y = 'SampleID')

grp = dplyr::group_by(tb, Genotype, Tissue, gid)
tb2 = dplyr::summarise(grp, nref = sum(nref), nalt = sum(nalt), ncft = sum(ncft))
tb3 = tb2[tb2$nref + tb2$nalt >= 30,]

fo = file.path(dirw, "20.ase.tissue.tsv")
#write.table(tb3, fo, sep = "\t", row.names = F, col.names = T, quote = F)

tb4 = ddply(tb3, .(Genotype, Tissue), summarise, ngene = length(gid))
tb5 = merge(tb3, tb4, by = c("Genotype", "Tissue"))
tb5$Genotype[tb5$Genotype == "B73xMo17"] = "BxM"
tb5$Genotype[tb5$Genotype == "Mo17xB73"] = "MxB"
tb6 = within(tb5, {
    lab = sprintf("%s %s [%d]", Tissue, Genotype, ngene)
    rm(Genotype, Tissue, ngene)
})

p1 = ggplot(tb6) +
  geom_histogram(aes(x = nref/(nref+nalt)), bins = 50) +
  scale_x_continuous(name = 'Proportion B73-specific Reads', limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous(name = "# Genes") +
  #scale_color_manual(name = "", values = cols) +
  facet_wrap( ~ lab, nrow = 4) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(legend.position = c(0.4, 0.8), legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = file.path(diro, "21.ase.tissue.pdf")
ggsave(p1, filename = fp, width = 12, height = 8)


# determine cis- or trans- regulatory pattern
fi = file.path(dirw, "20.ase.tissue.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)
ti2 = ti[ti$Genotype == 'B73xMo17', -1]

ff = file.path(diro, "../35.long.tsv")
tf = read.table(ff, sep = "\t", header = T, as.is = T)
tf2 = spread(tf[tf$Genotype %in% c("B73", "Mo17"),-4], Genotype, fpkm)

ti3 = merge(ti2, tf2, by = c("Tissue", "gid"))

# no way to tell cis or trans if both parents have 0 expression
idxs = which(ti3$B73 >= 1 | ti3$Mo17 >= 1)
ti4 = within(ti3[idxs,], { 
	c.pctB = B73 / (B73+Mo17)
	t.pctB = 0.5
})

c.pvals = apply(ti4, 1, myfunc <- function(x) 
	binom.test(as.integer(x[3]), as.integer(x[3])+as.integer(x[4]), as.numeric(x[9]), alternative = "two.sided")$p.value)
t.pvals = apply(ti4, 1, myfunc <- function(x) 
	binom.test(as.integer(x[3]), as.integer(x[3])+as.integer(x[4]), as.numeric(x[8]), alternative = "two.sided")$p.value)

ti5 = cbind(ti4, c.pval = c.pvals, t.pval = t.pvals)
#ti5 = cbind(ti4, c.pval = p.adjust(c.pvals, method = 'BH'), t.pval = p.adjust(t.pvals, method = 'BH'))

ti6 = data.frame(ti5, Reg = 'none', stringsAsFactors = F)
ti6$Reg[ti6$c.pval >= 0.05 & ti6$t.pval < 0.05] = 'cis'
ti6$Reg[ti6$c.pval < 0.05 & ti6$t.pval >= 0.05] = 'trans'
ti6$Reg[ti6$c.pval >= 0.05 & ti6$t.pval >= 0.05] = 'cis+trans'
ti6 = within(ti6, {
	reads = nref + nalt
	ptotal = B73 + Mo17
	p.propB = B73 / (B73+Mo17)
	h.propB = nref / (nref + nalt)
})
ti6$reads[ti6$reads > 1000] = 1000
ti6$ptotal[ti6$ptotal > 1000] = 1000
ti6$Reg[ti6$Reg == 'none' & ti6$p.propB < ti6$h.propB & ti6$h.propB > 0.5] = 'B73-biased'
ti6$Reg[ti6$Reg == 'none' & ti6$p.propB > ti6$h.propB & ti6$h.propB < 0.5] = 'Mo17-biased'

reg_levels = c("cis", "trans", "cis+trans", "B73-biased", "Mo17-biased", "none")
cols = c(brewer.pal(5, "Set1"), 'grey')
ti6$Reg = factor(ti6$Reg, levels = reg_levels)

# add DE
ff = file.path(dirw, '../43.deseq/01.de.tsv')
tf = read.table(ff, header = T, sep = "\t", as.is = T)
tf = tf[tf$comp == 'B73 vs Mo17', -2]
colnames(tf)[1] = "Tissue"

ti7 = merge(ti6, tf, by = c("Tissue", "gid"), all.x = T)
ti7$is.de[is.na(ti7$is.de)] = 'non-DE'

fo = file.path(dirw, "30.ase.tsv")
#write.table(ti7[,c('Tissue','gid','p.propB','h.propB','Reg','reads')], fo, sep = "\t", row.names = F, col.names = T, quote = F)


#### case study - ear tissue
ti8 = ti7[ti6$Tissue == 'ear_v14',]

# DE distribution
p1 = ggplot(ti8) +
  geom_point(aes(x = p.propB, y = h.propB, color = is.de)) +
  geom_hline(yintercept = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(name = 'Parental B73 Proportion', expand = c(0.01, 0)) +
  scale_y_continuous(name = 'F1 B73 Proportion', expand = c(0.01, 0)) +
  scale_color_brewer(name = '', palette = "Set1") +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(legend.position = c(0.4, 0.8), legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = file.path(dirw, "31.DE.pdf")
ggsave(p1, filename = fp, width = 8, height = 7)

# DE genes
p1 = ggplot(ti8[ti8$is.de != 'non-DE',]) +
  geom_point(aes(x = p.propB, y = h.propB, color = Reg)) +
  geom_hline(yintercept = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(name = 'Parental B73 Proportion', expand = c(0.01, 0)) +
  scale_y_continuous(name = 'F1 B73 Proportion', expand = c(0.01, 0)) +
  scale_color_manual(name = '', values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(legend.position = c(0.4, 0.8), legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = file.path(dirw, "32a.pdf")
ggsave(p1, filename = fp, width = 8, height = 7)

p1 = ggplot(ti8[ti8$is.de != 'non-DE',]) +
  geom_point(aes(x = p.propB, y = h.propB, color = Reg, size = reads)) +
  geom_hline(yintercept = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(name = 'Parental B73 Proportion', expand = c(0.01, 0)) +
  scale_y_continuous(name = 'F1 B73 Proportion', expand = c(0.01, 0)) +
  scale_color_manual(name = '', values = cols) +
  scale_size() +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(legend.position = c(0.4, 0.8), legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = file.path(dirw, "32b.pdf")
ggsave(p1, filename = fp, width = 8, height = 7)


# non-DE genes
p1 = ggplot(ti8[ti8$is.de == 'non-DE',]) +
  geom_point(aes(x = p.propB, y = h.propB, color = Reg)) +
  geom_hline(yintercept = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(name = 'Parental B73 Proportion', limits = c(0,1), expand = c(0.01, 0)) +
  scale_y_continuous(name = 'F1 B73 Proportion', limits = c(0,1), expand = c(0.01, 0)) +
  scale_color_manual(name = '', values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(legend.position = c(0.4, 0.8), legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = file.path(dirw, "33.pdf")
ggsave(p1, filename = fp, width = 8, height = 7)


## tissue-level summary
ti7s = ti7[ti7$is.de != 'non-DE',]
p1 = ggplot(ti7s) +
  geom_histogram(aes(x = Reg, fill = Reg), stat = "count") +
  #scale_x_d(name = 'Parental B73 Proportion', expand = c(0.01, 0)) +
  #scale_y_continuous(name = 'F1 B73 Proportion', expand = c(0.01, 0)) +
  scale_fill_manual(name = '', values = cols) +
  facet_wrap( ~ Tissue, nrow = 4) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 30, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = file.path(dirw, "35.cnt.pdf")
ggsave(p1, filename = fp, width = 8, height = 6)

