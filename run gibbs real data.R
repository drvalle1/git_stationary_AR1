rm(list=ls(all=TRUE))
library('mvtnorm')
library('Matrix')
set.seed(1)

setwd('U:/independent studies/kaplan amazon hydrology/gibbs stationarity/git_stationary_AR1')
source('gibbs functions.R')
source('gibbs_main.R')

setwd('U:/independent studies/kaplan amazon hydrology/hydrology data')
tmp=read.csv('final edited.csv',as.is=T)

#eliminate lagged variables
ind=grep('_lag',colnames(tmp)); colnames(tmp)[ind]
tmp1=tmp[,-ind]

#remove water level data after impact (i.e., after 1984)
ind=grep('w',colnames(tmp1)); colnames(tmp1)[ind]
cond=tmp1$ano >= 1984
tmp1[cond,ind]=NA
image(data.matrix(tmp1))

ngibbs=10000
a1.sigma2=b1.sigma2=0.1
ind=which(colnames(tmp1)%in%c('ano','mes'))
nomes=colnames(tmp1)[-ind]; nomes

for (i in 1:length(nomes)){
  dat=data.frame(y=tmp1[,nomes[i]],mes=tmp1$mes,ano=tmp1$ano)
  model=gibbs_main(dat=dat,ngibbs=ngibbs,a1.sigma2=a1.sigma2,b1.sigma2=b1.sigma2)
  
  #output results
  setwd('U:/independent studies/kaplan amazon hydrology/gibbs stationarity/results')
  nomes1=paste(c('betas_','outros_'),nomes[i],'.csv',sep='')
  seq1=seq(from=ngibbs*0.5,to=ngibbs,length.out=1000)
  write.csv(model$betas[seq1,],nomes1[1],row.names=F)
  write.csv(model$outros[seq1,],nomes1[2],row.names=F)
}



