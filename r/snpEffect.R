des = read.table(file.path(DIR_Out, 'vnt2_effect.txt'), header=TRUE);
des = cbind(des, 1);
cntCols = length(colnames(des));
colnames(des)[cntCols] = 'cnt';

cnames = colnames(des);
cmpIdxStart = 5;
cmpIdxEnd = 19;
rst = data.frame(0, 0, 0, 0, 0);
colnames(rst) = c('cmp', '-', 'ns', 's', 'ratio');
for (colIdx in cmpIdxStart:cmpIdxEnd) {
  cname = cnames[colIdx];
  w = aggregate(des$cnt, by=list(des[,colIdx]), FUN=sum);
  rowIdx = colIdx-cmpIdxStart+1;
  rst[rowIdx, ] = c(cname, as.numeric(w[,2]), 0);
}
rst$ns = as.numeric(rst$ns);
rst$s  = as.numeric(rst$s);
rst$ratio = rst$ns / (rst$ns+rst$s);
