# source("http://bioconductor.org/biocLite.R")
# biocLite("HTSFilter")
# biocLite("edgeR")

rm(list=ls())
setwd('/home/robin/Bureau/1507-Angers/TDcluster')
library(HTSCluster)
op = options()

#################################################################
# Data sets
#################################################################
# Drosophila, Subset data: embryos only
DataName = 'modencode_fly_pooled'
load(paste('../Data/', DataName,'.RData', sep=''))
dat <- exprs(modencodefly.eset.pooled)
conds <- pData(modencodefly.eset.pooled)
Cond <- as.vector(conds$stage[1:12])
X <- dat[,1:12]
X <- X[which(rowSums(X)>0),]

n = nrow(X); p = ncol(X)
hist(X, breaks=sqrt(prod(dim(X))))

#################################################################
# Clustering
#################################################################
# Poisson mixture
Kmax = 15
# PMC = PoisMixClusWrapper(y=X, conds=Cond, 
#                    gmin=1, gmax=Kmax, lib.type = "TC") 
load(paste(DataName, '_PMC.Rdata', sep=''))
# plot(PMC)
options(op)
Res = PMC$DDSE.results
K = length(Res$pi)
lambda = Res$lambda
pi = Res$lambda
par(pch=20)
plot(0, 0, xlim=c(1, p), ylim=c(min(lambda), max(lambda)), col=0, xlab='', ylab='')
for (k in 1:K){lines(lambda[, k], col=k, lwd=2)}

