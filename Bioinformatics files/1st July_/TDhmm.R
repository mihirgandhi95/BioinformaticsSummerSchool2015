rm(list=ls())

setwd('/home/robin/Bureau/1507-Angers/TDcluster')
source('/home/robin/PgmR/Mixture/FunctionsMclust.R')
library(RHmm)
library(mclust)
par(pch=20)

#################################################################
# Data sets
#################################################################
# Popova data
DataName = 'R_BLC_B1_T45-Chr16'
load(paste('../Data/', DataName, '.RData', sep=''))
Y = A+B

# Description
n = length(Y)
plot(Y)

#################################################################
# Mixture
#################################################################
hist(Y, breaks=sqrt(n))
GM = Mclust(Y)
# plot(GM)
F_PlotMclust1(Y, GM)

#################################################################
# HMM 
#################################################################
Kmax = 5
BIC = rep(NaN, Kmax); LogL = BIC
HMM = list()
LogL[1] = sum(dnorm(Y, mean=mean(Y), sd=sd(Y), log=T))
BIC[1] = -2*LogL[1]+2*log(n)
for (K in 2:Kmax){
   cat(K, '')
#    control = list(iter=1000)
   HMM[[K]] = HMMFit(Y, nStates=K)
   BIC[K] = HMM[[K]]$BIC
   LogL[K] = HMM[[K]]$LLH
}
#save(HMM, file=paste(DataName, '_HMM.Rdata', dep=''))
# Model selection
plot(BIC)
K = which.min(BIC)
HMMopt = HMM[[K]]
# Classification
V = viterbi(HMMopt, Y)
plot(Y, col=V$states)
FB = forwardBackward(HMMopt, Y)
Z = apply(FB$Gamma, 1, which.max)
table(V$states, Z)
