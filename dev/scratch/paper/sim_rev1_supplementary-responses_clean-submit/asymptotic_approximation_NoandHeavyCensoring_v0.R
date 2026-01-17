# Down 2 levels relative to project directory
source("../../R/source_forestsearch_v0.R")
source_fs_functions(file_loc="../../R/")

hr1 <- 1.25
hr2 <- 1.0

library(cubature)
library(kableExtra)
library(knitr)
library(ggplot2)

# No censoring
fileout <- c("results-input/pAnyH_approx_nsg_nocensoring.Rdata")
pC <- 0.0
hrH.plims<-seq(0.5,3.0,length=60)
hrH.plims <- sort(c(hrH.plims,0.5,0.7,0.75,0.80))

pAnyH.approx.100<-pAnyH.approx.80<-pAnyH.approx.60<-rep(NA,length(hrH.plims))
for(hh in 1:length(hrH.plims)){
  hrH.plim <- hrH.plims[hh]
  pAnyH.approx.60[hh] <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrH.plim,
  pC=pC,n.sg=60,k1=log(hr1),k2=log(hr2))$integral
  pAnyH.approx.80[hh] <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrH.plim,
  pC=pC,n.sg=80,k1=log(hr1),k2=log(hr2))$integral
  pAnyH.approx.100[hh] <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrH.plim,
  pC=pC,n.sg=100,k1=log(hr1),k2=log(hr2))$integral
}

save(hrH.plims,pAnyH.approx.60,pAnyH.approx.80,pAnyH.approx.100,file=fileout)


plot(hrH.plims,pAnyH.approx.60,xlab="Hazard ratio for subgroup H",ylab="Probability of any H found",type="l",lty=1,lwd=2.0,ylim=c(0,1))
abline(h=0.10,lwd=0.5,col="red")
lines(hrH.plims,pAnyH.approx.80,type="l",lty=1,lwd=2,col="grey")
lines(hrH.plims,pAnyH.approx.100,type="l",lty=1,lwd=2,col="blue")


# Heavy censoring
fileout <- c("results-input/pAnyH_approx_nsg_heavycensoring.Rdata")
pC <- 0.80
hrH.plims<-seq(0.5,3.0,length=60)
hrH.plims <- sort(c(hrH.plims,0.5,0.7,0.75,0.80))

pAnyH.approx.100<-pAnyH.approx.80<-pAnyH.approx.60<-rep(NA,length(hrH.plims))
for(hh in 1:length(hrH.plims)){
  hrH.plim <- hrH.plims[hh]
  pAnyH.approx.60[hh] <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrH.plim,
                                        pC=pC,n.sg=60,k1=log(hr1),k2=log(hr2))$integral
  pAnyH.approx.80[hh] <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrH.plim,
                                        pC=pC,n.sg=80,k1=log(hr1),k2=log(hr2))$integral
  pAnyH.approx.100[hh] <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrH.plim,
                                         pC=pC,n.sg=100,k1=log(hr1),k2=log(hr2))$integral
}

save(hrH.plims,pAnyH.approx.60,pAnyH.approx.80,pAnyH.approx.100,file=fileout)

plot(hrH.plims,pAnyH.approx.60,xlab="Hazard ratio for subgroup H",ylab="Probability of any H found",type="l",lty=1,lwd=2.0,ylim=c(0,1))
abline(h=0.10,lwd=0.5,col="red")
lines(hrH.plims,pAnyH.approx.80,type="l",lty=1,lwd=2,col="grey")
lines(hrH.plims,pAnyH.approx.100,type="l",lty=1,lwd=2,col="blue")




