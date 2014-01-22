source(file.path("/Volumes/pzhou/Scripts/r", "my.r"));

chrName = "MtChr2";
fWindow = file.path(DIR_Stat, chrName, "window.txt");
win = read.table(fWindow, header=FALSE);
win = cbind(win, 0, 0, 0, 0);
colnames(win) = c("cnt", "snpStart", "snpStop", "chrStart", "chrStop", "cntSnp", "length", "posMiddle", "rhoMean", "rhoUp", "rhoDown");
win$posMiddle = ((win$chrStop + win$chrStart) / 2 ) / 1000;
for (i in 1:400) {
  fRate = file.path(DIR_Stat, chrName, sprintf("rates_%03d.txt", i));
  fLoc = file.path(DIR_Misc, chrName, sprintf("in_LDhat_%03d_loc.txt", i));
  summary <- summarise.interval(rates.file=fRate, locs.file=fLoc);
  win$rhoMean[i] = mean(summary[1,1]);
  win$rhoUp[i] = mean(summary[1,2]);
  win$rhoDown[i] = mean(summary[1,4]);
}
write.table(win, file = file.path(DIR_Stat,paste("rho_",chrName,".txt",sep="")), sep="\t", quote=FALSE);
p <- ggplot(win) +
  geom_line(size=0.3, aes(posMiddle, rhoMean), colour="red") +
  geom_line(size=0.2, aes(posMiddle, rhoUp), colour="blue") +
  geom_line(size=0.2, aes(posMiddle, rhoDown), colour="green") +
  xlab("Location on chromosome (/100kb)") + ylab("rho") +
  opts(title=chrName) +
  opts(axis.text.x = theme_text(size = 7, colour = "grey50"));
ggsave(p, filename = file.path(DIR_Stat,paste("rho_",chrName,".png",sep="")),
    width=8, height=5);