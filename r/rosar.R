require(plyr)
require(dplyr)
require(ggplot2)
source('Location.R')

dirw = file.path(Sys.getenv("rosar"), "test")

fi = file.path(dirw, "21.tbl")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
