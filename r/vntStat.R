source("/Volumes/pzhou/Scripts/r/my.r");

d <- readIn("geneStat.txt");
t = agg1(d);
#plotT(t);
g <- agg2(d);
#plotG(g);
#plotG_LR(g);
a <- agg3(d);

subGeneNames = c("Hydrolase(299)", "Late_nodulin(236)", "Peptidase(305)", 
  "Protein_kinase(569)");
subAccs = c("HM002", "HM015", "HM017", "HM101");
#plotA(a, subGeneNames, subAccs);
#plotA(a);

chrLst = 1:8;
subAccs = c("HM002", "HM015", "HM017", "HM029", "HM101");
for (i in chrLst) {
  chrName = paste("MtChr",i,sep="");
  fName = file.path(DIR_Stat, "chrStat_10k", paste("stat_",chrName,".txt",sep=""));
#  plotChr(fName, subAccs);
}

fName = "geneDesc_0102_snp.txt";
subGeneNames = c("Hydrolase(299)", "Late_nodulin(236)", "Peptidase(305)",
  "Protein_kinase(569)");
subAccs = c("HM002", "HM017", "HM015");
option = 1;
#plotSNPEffect(fName, option, subGeneNames, subAccs);
#plotSNPEffect(fName, option);

#plotRsq("rsq.out");