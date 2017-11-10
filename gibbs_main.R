gibbs_main=function(dat,ngibbs,a1.sigma2,b1.sigma2){
  #useful
  ntime=nrow(dat)
  aux=outer(1:ntime,1:ntime,'-')
  aux1=abs(aux)
  numeros=matrix(1:(ntime^2),ntime,ntime)
  lo.diag=diag(numeros)+1; lo.diag=lo.diag[-length(lo.diag)]
  up.diag=diag(numeros)-1; up.diag=up.diag[-1]
  
  #add seasonal effects
  xmat=numeric()
  for (i in 1:12){
    cond=dat$mes==i
    tmp=ifelse(cond,1,0)
    xmat=cbind(xmat,tmp)
  }
  colnames(xmat)=paste('mes',1:12,sep='')
  xmat=cbind(xmat,dat$tempo,dat$tempo^2)
  
  #remove NA's
  y=dat$y
  ind.nas=which(is.na(y))
  y1=y[-ind.nas]
  xmat1=Matrix(xmat[-ind.nas,])
  aux2=aux1[-ind.nas,-ind.nas]
  
  #initial values
  betas=rep(0,ncol(xmat))
  sigma2=1
  rho=0.7
  Sigma=rho^aux2
  invSigma=solve(Sigma)
  
  #gibbs stuff
  vec.betas=matrix(NA,ngibbs,ncol(xmat))
  vec.outros=matrix(NA,ngibbs,2)
  colnames(vec.outros)=c('sigma2','rho')
  
  param=list(betas=betas,sigma2=sigma2,rho=rho,Sigma=Sigma,
             invSigma=Matrix(invSigma))
  jump1=0.05
  accept1=0
  
  for (i in 1:ngibbs){
    print(i)
    param$betas=update.betas(param,xmat1,y1)#betas.true#
    param$sigma2=update.sigma2(param,y1,a1.sigma2,b1.sigma2,xmat1)
    
    tmp=update.rho(param,jump1,y1,xmat1,aux2)
    param$rho=tmp$rho #rho.true#
    param$Sigma=tmp$Sigma #rho.true^aux2#
    param$invSigma=tmp$invSigma #inv.ar1(rho.true)[-ind.nas,-ind.nas] #
    accept1=accept1+tmp$accept
    
    #MH adaptation
    if (i%%50==0 & i<500){
      z=accept1/50
      print(z); print(jump1)
      if (z>0.3 & jump1<3)     jump1=jump1*2
      if (z<0.1 & jump1>0.001)  jump1=jump1*0.5
      accept1=0  
    }
    
    vec.betas[i,]=param$betas
    vec.outros[i,]=c(param$sigma2,param$rho)
  }
  
  list(betas=vec.betas,outros=vec.outros)  
}

