require(rtracklayer)
source("comp.fun.R")

tname = "hm101"
qname = "hm056"
rname = "hm340"

t = read_genome_stat(tname)
q = read_genome_stat(qname)
r = read_genome_stat(rname)

vq = read_var_stat(qname)
vr = read_var_stat(rname)

cq = read_comp_stat(qname, tname)
cr = read_comp_stat(rname, tname)

# basic assembly statistics
source(paste("http://faculty.ucr.edu/~tgirke/",
    "Documents/R_BioCond/My_R_Scripts/contigStats.R", sep=''))
library(Biostrings)
assembly <- readDNAStringSet(
    sprintf("%s/11_genome.fa", q$dir, q$name), "fasta")
N <- list(acc = width(assembly))
reflength <- sapply(N, sum)
stats <- contigStats(N = N, reflength = reflength, style = "data")
stats[["Contig_Stats"]]

# plot scaffold size distribution
tmp = cut(t_len$length / 1000, breaks = c(0,1,5,10,50,100,500,1000,5000))
p = ggplot(data.frame(size = tmp)) +
  geom_bar(aes(x = factor(size)), width = 0.7) + 
  scale_x_discrete(name = "Scaffold Size (kb)") + 
  scale_y_continuous(name = "") +
  theme(axis.text.x = element_text(angle = 15, size = 8))
ggsave(file.path(q$dir, "figs/01_scaffold_size.png"), p, width=5, height=4)

# plot global pairwise comparison
p <- ggplot(data = q$tw) +
  geom_rect(mapping = aes(xmin = hBeg / 1000000, xmax = hEnd/1000000, 
    ymin = 0, ymax = 1, fill = hId)) +  
  layer(data = tg, geom = 'rect', mapping = 
    aes(xmin = hBeg / 1000000, xmax = hEnd / 1000000, ymin = -1, ymax = 0), 
        geom_params = list()) +
  scale_x_continuous(name = 'Chr Position (Mbp)', expand = c(0.01, 0)) +
  scale_y_continuous(name  ='', expand = c(0.04, 0)) +
  facet_grid(hId ~ .) + 
  theme(legend.position = 'right', legend.title = element_blank()) +
  theme(axis.text.x = element_text(siz e= 8, angle = 0)) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank())
ggsave(p, filename = file.path(q$dir, "figs/03_coverage.png"), 
    width = 7, height = 5)
