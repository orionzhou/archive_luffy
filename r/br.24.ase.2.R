require(plyr)
require(dplyr)
require(tidyr)
require(ggplot2)
require(RColorBrewer)

dirw = file.path(Sys.getenv("misc2"), "briggs")
dirw = '/home/springer/zhoux379/scratch/briggs'
diro = file.path(Sys.getenv("misc2"), "briggs", "46.ase")

fl = file.path(dirw, '00.1.read.correct.tsv')
tl = read.table(fl, header = T, sep = "\t", as.is = T)[,1:5]
tl = within(tl, {label = sprintf("%s: %s %s rep%d", SampleID, Tissue, Genotype, Treatment)})


## Indel mapping bias
sids = c("BR006", "BR054", "BR151")
opts = c("snp+indel", "snp")

for (sid in sids) {
for (opt in opts) {


fi = sprintf("%s/25.ase/%s.5.%s.bed", dirw, sid, opt)
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

## calc gene stats
for (i in 101:nrow(tl)) {
for (i in 61:76) {
	sid = tl$SampleID[i]
	f1 = sprintf("%s/25.ase/%s.tsv", dirw, sid)
	f2 = sprintf("%s/25.ase/%s.bed", dirw, sid)
	if(! file.exists(f1) | !file.exists(f2)) {
		cat(sid, 'skipped\n')
		next
	}
	fo = sprintf("%s/01.all/%s.tsv", diro, sid)
	cmd = sprintf("bed.ase.sum.py %s %s %s", f1, f2, fo)
	system(cmd)
	cat(sid, "\n")
}

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
	ti = ti[ti$nref + ti$nalt >= 30, ]
	
	geno = ifelse(geno == "Mo17xB73", "MxB", "BxM")
	lab = sprintf("%s %s [%d]", tiss, geno, nrow(ti))
	
	ti = cbind(sid = sid, lab = lab, ti)
	to = rbind(to, ti)
	cat(sid, "\n")
}
fo = file.path(diro, "10.ase.gene.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)


# plot ASE summary
fi = file.path(diro, "10.ase.gene.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

p1 = ggplot(ti) +
  geom_histogram(aes(x = nref/(nref+nalt)), bins = 50) +
  scale_x_continuous(name = 'Proportion B73-specific Reads', limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous(name = "# Genes") +
  #scale_color_manual(name = "", values = cols) +
  facet_wrap( ~ lab, nrow = 8) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  #theme(legend.position = c(0.4, 0.8), legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = file.path(diro, "11.ase.gene.pdf")
ggsave(p1, filename = fp, width = 13, height = 13)
