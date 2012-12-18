pre = 3
opt = 'acc26'

pre = 2 
opt = 'acc56'

pre = 0 
opt = 'acc84'

accs = get_mt_ids(opt)
dirI = sprintf("%s/mt_35/31_phylogeny/%02d/08_stat", DIR_Repo, pre)
cutoff_missing = length(accs) * 0.3
intervals = seq(0,0.5,0.05)

chr = 5
fi = file.path(dirI, paste("chr", chr, ".tbl", sep=''))
s01 = read.table(fi, header=T, sep="\t", as.is=T)
s02 = cbind(s01, freq_der = s01$n_der / (s01$n_anc+s01$n_der))
s03 = s02[s02$n_states == 2 & s02$n_N < cutoff_missing,]
t01 = table(cut(s03$freq_der, breaks=intervals))

df = data.frame(bin=names(t01), count=as.numeric(t01))
p = ggplot(df) +
	geom_bar(aes(x=bin, y=count, fill='all', width=0.6), stat='identity', position='dodge') + 
	scale_fill_brewer(palette='Set3') +
	scale_x_discrete(name="Minor Alelle Frequency") +
	scale_y_continuous(name="Number SNPs", formatter="comma") +
	opts(title=paste("MAF distribution of chr", chr, " SNPs", sep=""), axis.text.x = theme_text(angle=45, size=8))
fo = file.path(dirI, paste("chr", chr, "_sfs.png", sep=""))
ggsave(fo, p, width=5, height=4)



