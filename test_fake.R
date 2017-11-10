# rm(list=ls(all=TRUE))
library('mvtnorm')
library('Matrix')
set.seed(1)

setwd('U:/independent studies/kaplan amazon hydrology/gibbs stationarity/git_stationary_AR1')
source('gibbs functions.R')
source('gibbs_main.R')
dat=read.csv('fake data.csv',as.is=T)
dat1=data.frame(y=dat[,1],mes=dat$mes,tempo=dat$tempo)

ngibbs=2000
a1.sigma2=b1.sigma2=0.1
model=gibbs_main(dat=dat1,ngibbs=ngibbs,a1.sigma2=a1.sigma2,b1.sigma2=b1.sigma2)

#check to see if it works
check=function(true,estim){
  rango=range(c(true,estim))
  plot(true,estim,xlim=rango,ylim=rango)
  lines(rango,rango)
}
seq1=(ngibbs/2):ngibbs
betas.estim=colMeans(model$betas[seq1,])
check(betas.true,betas.estim)

plot(model$outros[seq1,'sigma2'],type='l')
abline(h=sigma2.true,col='blue')

plot(model$outros[seq1,'rho'],type='l')
abline(h=rho.true,col='blue')