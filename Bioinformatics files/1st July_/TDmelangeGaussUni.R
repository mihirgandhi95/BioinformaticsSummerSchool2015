# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")

rm(list=ls())
setwd('/home/robin/Bureau/1507-Angers/TDcluster')
library(mclust)
source('/home/robin/PgmR/Mixture/FunctionsMclust.R')

#################################################################
# Data set
#################################################################
# Nematode data 
data = read.table('../Data/datasetNematode.dat', h=T)
Gene = data[, 1]
Y = t(log2(data[, -1]))

# Descriptive stats
n = nrow(Y); p = ncol(Y)
hist(apply(Y, 1, sd), breaks=sqrt(n))
hist(apply(Y, 1, mean), breaks=sqrt(n))

#################################################################
# Clustering for 1 var
#################################################################
c = 26
y = Y[, c]
GM = Mclust(y)
plot(GM)
F_PlotMclust1(y, GM)

