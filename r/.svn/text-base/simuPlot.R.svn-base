require(ggplot2);
if(.Platform$OS.type=='windows') {
	dir <- file.path(Sys.getenv("HOMEPATH"), "Scripts", "out");
} else {
	dir <- file.path(Sys.getenv("HOME"), "Scripts", "out");
}
dataf <- read.table(file.path(dir, "out4R.dat"), header=TRUE);
SubPop = paste(c(""),t(dataf["pop"]+1),sep="");
#png(filename = file.path(dir, "rg.png"));
p <- ggplot(dataf, aes(gen, af)) + 
	xlab("Generation") + ylab("Allele Frequency") +
	ylim(0,1) +
	aes(shape=SubPop) + 
	geom_point(aes(colour=SubPop)) + 
	geom_line(aes(colour=SubPop)) + 
	opts(title = "Random Drift");
ggsave(p, filename = file.path(dir,"rg.png"), width=10, height=5);
#dev.off();
	
#qplot(gen, af, data = dataf, colour = pop, 
# shape = paste(c("Subpop"),t(dataf["pop"]),sep=""), group = pop, geom="point", 
#	xlab="Generation", ylab="Allele Frequency", main="Random Drift", 
#	legend = "Subpop",
#	ylim=c(0,1)
#);
