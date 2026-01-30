
dens_threshold_both<-function(x,theta,pC=0.3,n.sg,k1,k2){
  d.sg<-n.sg*(1-pC)
  mu_split<-log(theta)
  sig2_split<-8/d.sg
  sig_split<-sqrt(sig2_split)
  dens_1<-dnorm(x[1],mean=mu_split,sd=sig_split)
  dens_2<-dnorm(x[2],mean=mu_split,sd=sig_split)
  ans<-c((x[1]+x[2]) >= 2*k1)*c(x[1]>=k2)*c(x[2]>=k2)*dens_1*dens_2
  return(ans)
}


hr1 <- 1.25
hr2 <- 1.0

library(cubature)

# Mild censoring

n.sg <- 60
pC <- 0.20


hrH.plims<-seq(0.5,3.0,length=60)
hrH.plims <- sort(c(hrH.plims,0.5,0.7,0.75,0.80))


pAnyH.approx <-rep(NA,length(hrH.plims))
for(hh in 1:length(hrH.plims)){
  hrH.plim <- hrH.plims[hh]
  pAnyH.approx[hh] <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrH.plim,
  pC=pC, n.sg = n.sg, k1=log(hr1),k2=log(hr2))$integral
}

plot(hrH.plims,pAnyH.approx,xlab="Hazard ratio for subgroup H",ylab="Probability of any H found",type="l",lty=1,lwd=2.0,ylim=c(0,1))
abline(h=0.10,lwd=0.5,col="red")
