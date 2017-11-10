tnorm <- function(lo,hi,mu,sig){ 
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(1,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#----------------------------------------
fix.MH=function(lo,hi,old1,new1,jump){
  jold=pnorm(hi,mean=old1,sd=jump)-pnorm(lo,mean=old1,sd=jump)
  jnew=pnorm(hi,mean=new1,sd=jump)-pnorm(lo,mean=new1,sd=jump)
  log(jold)-log(jnew) #add this to pnew
}
#----------------------------------------------------------------------------------------------
acceptMH <- function(p0,p1,x0,x1){   #accept for M, M-H
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(1,0,1)
  if (z<a) return(x1)
  return(x0)
}
#---------------------------------
update.betas=function(param,xmat1,y1){
  k=1/param$sigma2
  prec=data.matrix(t(xmat1)%*%(k*param$invSigma)%*%xmat1)+diag(1,ncol(xmat1))
  var1=solve(prec)
  pmedia=t(xmat1)%*%(k*param$invSigma)%*%y1
  t(rmvnorm(1,mean=var1%*%pmedia,sigma=var1))
}
#---------------------------------------------------------
update.sigma2=function(param,y1,a1.sigma2,b1.sigma2,xmat1){
  a1=(length(y1)/2)+a1.sigma2
  err=y1-xmat1%*%param$betas
  pb=t(err)%*%param$invSigma%*%err
  b1=(1/2)*pb+b1.sigma2
  1/rgamma(1,a1,as.numeric(b1))
}
#---------------------------------------------------------
update.rho=function(param,jump,y1,xmat1,aux2){
  err=y1-xmat1%*%param$betas
  rho.old=param$rho
  p1.old=determinant(param$Sigma,logarithm = T)$modulus[[1]]
  p2.old=(1/param$sigma2)*t(err)%*%param$invSigma%*%err
  
  rho.new=tnorm(lo=0,hi=1,mu=rho.old,sig=jump)
  Sigma.new=rho.new^aux2
  invSigma.new=solve(Sigma.new)
  p1.new=determinant(Sigma.new,logarithm = T)$modulus[[1]]
  p2.new=(1/param$sigma2)*t(err)%*%invSigma.new%*%err
  ajuste=fix.MH(lo=0,hi=1,rho.old,rho.new,jump)
  
  k=acceptMH(-(1/2)*as.numeric(p1.old+p2.old),
             -(1/2)*as.numeric(p1.new+p2.new)+ajuste,
             rho.old,rho.new)
  
  if (k==rho.old) return(list(rho=rho.old,Sigma=param$Sigma,invSigma=param$invSigma,accept=0))
  if (k!=rho.old) return(list(rho=rho.new,Sigma=Sigma.new,invSigma=Matrix(invSigma.new),accept=1))
}
