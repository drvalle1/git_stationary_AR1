rm(list=ls(all=TRUE))
library('mvtnorm')
library('Matrix')
set.seed(1)

setwd('U:/independent studies/kaplan amazon hydrology/gibbs stationarity/git_stationary_AR1')
source('gibbs functions.R')
source('gibbs_main.R')

setwd('U:/independent studies/kaplan amazon hydrology/hydrology data')
tmp=read.csv('final edited.csv',as.is=T)

dat=read.csv('fake data.csv',as.is=T)
dat1=data.frame(y=dat[,1],mes=dat$mes,tempo=dat$tempo)
model=gibbs_main(dat1)

setwd('U:/independent studies/kaplan amazon hydrology/hydrology data/stationarity/results')
nomes=colnames(dat)[uuu]
nomes1=paste(c('betas_','outros_'),nomes,'.csv',sep='')
seq1=seq(from=ngibbs*0.5,to=ngibbs,length.out=1000)
write.csv(vec.betas[seq1,],nomes1[1],row.names=F)
write.csv(vec.outros[seq1,],nomes1[2],row.names=F)
}
