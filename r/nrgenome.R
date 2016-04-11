require(rtracklayer)
require(seqinr)
require(GenomicRanges)
source('Location.R')
source('comp.fun.R')

qnames = c(tname, qnames_ingroup)

qname = qnames[1]
dirw = file.path(Sys.getenv("misc2"), "nrgenome", qname)

##### create sliding window table
tlen = read.table(cfgs[[qname]]$size, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

tgap = read.table(cfgs[[qname]]$gap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

winsize = 100
x = tt$end
names(x) = tt$chr
gr = tileGenome(x, tilewidth = winsize, cut.last.tile.in.chrom = T)

tw = data.frame(chr = seqnames(gr), beg = start(gr), end = end(gr), 
  len = width(gr), stringsAsFactors = F)
  
bp_gap = intersect_basepair(gr, grp)
idxs = bp_gap / tw$len > 0.5 | tw$len < 50
cat(sum(tw$len[idxs]) / 1000000, "Mbp removed\n")

to = tw[!idxs, 1:3]
to$beg = to$beg - 1
fo = file.path(dirw, "01.win100.bed")
write.table(to, fo, sep = "\t", row.names = F, col.names = F, quote = F)

## seqret.pl -d $genome/HM101/11_genome.fas -b 01.bed -o 01.fas
## qsub.blat.pl -i 01.fas -o 11.blat -n 1 -p 2 -t $genome/HM101/11_genome.fas

