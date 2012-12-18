fi = file.path(DRI_Misc1, "cat/mt_35/03_gene_go.tbl")
c01 = read.table(fi, header=TRUE, sep="\t")
dirO = file.path(DIR_Misc1, "r", "cat")

c11 = data.frame(table(c01$go_id))
colnames(c11) = c("go_id", "count")
c12 = cbind(c11, ontology=NA, desc=NA)
for (i in 1:nrow(c12)) {
	go = GOTERM[[ as.character(c12$go_id[i]) ]]
	if( !is.null(go) ) {
		c12$ontology[i] = Ontology(go)
		c12$desc[i] = Term(go)
	}
}
c12 = c12[order(-c12$count),]
write.table(c12, file.path(dirO, "11_go.tbl"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

c13 = droplevels(subset(c12, count >= 50))
c13$go_id = factor(c13$go_id, levels=c13$go_id[order(c13$count)])
p = ggplot(c13) + 
	geom_bar(aes(x=go_id, y=count, fill=ontology), width=0.7, stat='identity') +	scale_fill_brewer(palette='Set1') +
	scale_x_discrete(name='', breaks=c13$go_id, labels=c13$desc) + 
	coord_flip()
p
ggsave(p, filename = file.path(DIR_R, 'cat/12.png'), width=8, height=8);
