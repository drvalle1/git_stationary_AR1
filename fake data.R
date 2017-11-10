rm(list=ls(all=TRUE))
library('mvtnorm')
set.seed(199)

setwd('U:/independent studies/kaplan amazon hydrology/gibbs stationarity')
dat=read.csv('data for stationarity.csv',as.is=T)

ntime=nrow(dat)
aux=outer(1:ntime,1:ntime,'-')
aux1=abs(aux)

#add seasonal effects
xmat=numeric()
for (i in 1:12){
  cond=dat$mes==i
  tmp=ifelse(cond,1,0)
  xmat=cbind(xmat,tmp)
}
colnames(xmat)=paste('mes',1:12,sep='')
xmat=cbind(xmat,dat$tempo,dat$tempo^2)

#parameter values
rho.true=rho=0.1
sigma2.true=sigma2=0.5
betas.true=betas=runif(ncol(xmat))

#generate y
media=xmat%*%betas
y=rmvnorm(1,mean=media,sigma=sigma2*(rho^aux1))
y1=t(y)

#put in some missing data
ind=sample(1:ntime,size=ntime*0.1)
y1[ind]=NA
dat[,1]=y1

setwd('U:/independent studies/kaplan amazon hydrology/gibbs stationarity/git_stationary_AR1')
write.csv(dat,'fake data.csv',row.names=F)

