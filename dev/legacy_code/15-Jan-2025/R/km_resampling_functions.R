
plot.band<-function(x,mean.value,lower,upper,show.axes=F,band=TRUE,ltype="l",lty=1,xlabel=NULL,ylabel=NULL,color="grey",ylim=c(min(lower,na.rm=TRUE),max(upper,na.rm=TRUE))){
plot(x[order(x)],mean.value[order(x)],type="n",axes=show.axes,xlab=xlabel,lty=lty,
ylab=ylabel,ylim=ylim)
if(band){
polygon(c(x[order(x)],rev(x[order(x)])),c(lower[order(x)],rev(upper[order(x)])),col=color,border=FALSE)
lines(x[order(x)],mean.value[order(x)],lty=lty,lwd=2.5,type=ltype)
}
}


plot.band.two<-function(x,curve1,curve2,lower,upper,show.axes=F,
ltype="l",lty=1,xlabel=NULL,ylabel=NULL,ylim=c(min(lower,na.rm=TRUE),max(upper,na.rm=TRUE))){
plot(x[order(x)],curve1[order(x)],type="n",axes=show.axes,xlab=xlabel,lty=lty,
ylab=ylabel,ylim=ylim)
polygon(c(x[order(x)],rev(x[order(x)])),
c(lower[order(x)],rev(upper[order(x)])),col="lightgrey",border=F)
lines(x[order(x)],curve1[order(x)],lty=lty,lwd=2.5,type=ltype)
lines(x[order(x)],curve2[order(x)],lty=lty,lwd=2.5,type=ltype)
}



N.Weighted.KM<-function(x,error,w=1){
sum(w*(error<=x))
}

R.Weighted.KM<-function(x,error,w=1){
sum(w*(error>=x))
}


NA.CHR.Weighted.RMST<-function(time,delta,w.n=rep(1,length(time)),w.d=rep(1,length(time)),
at.points=sort(time),se.type="greenwood",get.Stute=FALSE,tpoints.add=NULL){

if(!is.null(tpoints.add)) at.points<-sort(c(unique(c(at.points,tpoints.add))))

if(se.type!="greenwood" & se.type!="tsiatis") stop("Invalid se type -- greenwood or tsiatis allowed")
#is.sorted<-(all(time==sort(time)))
is.sorted<-!is.unsorted(time)
if(!is.sorted){
id<-order(time); time<-time[id]; delta<-delta[id]; w.n<-w.n[id]; w.d<-w.d[id]
}
risk<-unlist(lapply(as.list(at.points),R.Weighted.KM,error=time,w=w.d))
###########################################################################
### Adaptive H process correspoinding to N-A rep. via integral wrt M.G ####
Hmart.chf<-ifelse(risk>0,1/risk,0)
############################################################################
counting<-unlist(lapply(as.list(at.points),N.Weighted.KM,error=time,w=w.n*ifelse(delta==1,1,0)))
counting <- c(0, counting)
delta.counting <- diff(counting)
DN.Risk<-ifelse(risk>0,delta.counting/risk,0.0)
chf <- cumsum(DN.Risk)
var.chf<-cumsum(ifelse(risk>0,delta.counting/(risk^2),0.0))
S.KM <- cumprod(1-DN.Risk)

#if(any(S.KM<0)){
#cat("t,S(t),1-DN.risk,DN.risk","\n")
#print(cbind(at.points,S.KM,1-DN.Risk,DN.Risk))
#}

#S.KM[which(S.KM<0)]<-NA
S.KM[which(S.KM<0)]<-0.0

S.NA <- exp(-chf)
var.NA<-(S.NA^2)*var.chf
# Greenwood variance estimate
if(se.type=="greenwood"){
aa<-delta.counting
bb<-risk*(risk-delta.counting)
var.KM<-(S.KM^2)*cumsum(ifelse(risk>0,aa/bb,0.0))
se.KM<-sqrt(var.KM)
}
if(se.type=="tsiatis"){
var.KM<-(S.KM^2)*var.chf
se.KM<-sqrt(var.KM)
}
result<-list(time=time,at.points=at.points,S.NA=S.NA,S.KM=S.KM,chf=chf,se.chf=sqrt(var.chf),se.NA=sqrt(var.NA),
n.risk=risk,n.events=delta.counting,Hmart.chf=Hmart.chf,se.KM=se.KM)
return(result)
}




# Mean for general S (truncated at L)
# \int_{0}^{L}S(t)dt
# fatx=S, x=t
LS.int<-function(fatx,x,L=Inf){
f<-fatx[1:(length(fatx)-1)]
d.x<-diff(x)
int.f<-x[1]+sum((d.x*f)[which(x[-1]<=L)])
return(int.f)
}


mu.L<-function(L,S,t){
mu<-LS.int(fatx=S,x=t,L=L)
return(mu)
}

# Need to check!
LS.int.del<-function(fx,x,dx,L=Inf){
  int.f<-sum((dx*fx)[which(x<=L)])
  return(int.f)
}


#################################################
# Check LS.int() integration with normal density
# and the integrate(.) function
#################################################
#integrate(dnorm, 0, 200)
#integrate(dnorm, 0, 2000)
# check LS.int
#test1<-integrate(dnorm, 0, 2)
#x<-seq(0,2,by=0.0001)
#f1<-dnorm(x)
#test1.check<-LS.int(f1,x)
#print(c(test1$value,test1.check))

#test1<-integrate(dnorm, 0, 200)
#x<-seq(0,200,by=0.0001)
#f1<-dnorm(x)
#test1.check<-LS.int(f1,x)
#print(c(test1$value,test1.check))

#test1<-integrate(dnorm, 0, 2000)
#x<-seq(0,2000,by=0.0001)
#f1<-dnorm(x)
#test1.check<-LS.int(f1,x)
#print(c(test1$value,test1.check))

#test1<-integrate(dnorm, 0, 2000)
#x<-seq(0,3000,by=0.0001)
#f1<-dnorm(x)
#test1.check<-LS.int(f1,x,L=2000)
#print(c(test1$value,test1.check))

#test1<-integrate(dnorm, 0, 2)
#x<-seq(0,30,by=0.0001)
#f1<-dnorm(x)
#test1.check<-LS.int(f1,x,L=2)
#print(c(test1$value,test1.check))




surv.tau<-function(S,tpoints,tau){
tt<-c(0,tpoints)
ss<-c(1,S)
S.point<-approx(tt,ss,xout=c(tau),method='constant',f=0)$y
return(S.point)
}

surv.quantile<-function(S,tpoints,tau=0.5){
S.tau<-min(tpoints[S<=(1-tau)])
return(S.tau)
}

one<-function(x){
return(rep(1,length(x)))
}

zero<-function(x){
return(rep(0,length(x)))
}



##############################################################################
#Note:  tpoints.add<-{t1,t2,t3} allows for estimation at specific time points
# which may not be unique events
###############################################################################
KM.resample<-function(time,delta,w.n=rep(1,length(time)),w.d=rep(1,length(time)),
tpoints.add=NULL,draws=0,details=FALSE,at.points=NULL,G.draws=NULL){

if(is.null(at.points)) at.points<-sort(unique(time))

is.sorted<-!is.unsorted(time)
if(!is.sorted){
id<-order(time); time<-time[id]; delta<-delta[id]; w.n<-w.n[id]; w.d<-w.d[id]
}

if(!is.null(G.draws)) draws<-ncol(G.draws)

temp<-NA.CHR.Weighted.RMST(time=time,delta=(delta==1),w.n=w.n,w.d=w.d,tpoints.add=tpoints.add,at.points=at.points)
S.KM<-temp$S.KM
at.points<-temp$at.points
n<-length(time)
if(draws>0 | !is.null(G.draws)){
# Distribution of sqrt(n){\hat{S}(.)-S_{true}}
S.KM.star<-matrix(NA,nrow=length(at.points),ncol=draws)
if(is.null(G.draws)) G.draws<-matrix(rnorm(draws*n),ncol=draws)
for(dd in 1:draws){
g<-G.draws[,dd]
risk<-unlist(lapply(as.list(at.points),R.Weighted.KM,error=time,w=w.d))
Wn.g<-w.n*g
counting<-unlist(lapply(as.list(at.points),N.Weighted.KM,error=time[delta==1],w=Wn.g[delta==1]))
counting <- c(0, counting)
delta.counting <- diff(counting)
S.KM.star[,dd]<--S.KM*cumsum(ifelse(risk>0,delta.counting/risk,0.0))
if(details){
if(round(dd/draws,digits=3)==0.10) cat("10% resampling done","\n")
if(round(dd/draws,digits=3)==0.25) cat("25% resampling done","\n")
if(round(dd/draws,digits=3)==0.75) cat("75% resampling done","\n")
if(round(dd/draws,digits=3)==0.90) cat("90% resampling done","\n")
}
}
}
if(draws>0) out<-list(S.KM=S.KM,at.points=at.points,S.KM.star=S.KM.star)
if(draws==0) out<-data.frame(S.KM,at.points)
return(out)
}



# Note: In most code time indicates the observed survival times; delta event indicator, and z=group
# At some point, want to align!
# Especially as in other code time denotes the risk set! (below, R= risk set)

KM.MeanTrunc<-function(time,delta,z,L=NULL,draws=0,G1.draws=NULL,G0.draws=NULL,get.band=FALSE,taus=c(NULL,NULL),tau.seq=0.25,draws.band=draws,details=FALSE,plotband=get.band,L.type="events"){
###############################################################################
# Optional is the calculation of simultaneous bands (get.Band=TRUE)
# sup_{tau1,tau2}{(\hat{S1}-\hat{S0})-(S1_{true}-S0_{true})}/\hat{sig2}
###############################################################################

if(!is.null(G0.draws) & !is.null(G1.draws)) draws<-ncol(G0.draws)

if(is.null(L)){
################
# Control group
################
to.get<-which(z==0)
time0<-time[to.get]
delta0<-delta[to.get]
#####################
# Experimental group
####################
to.get<-which(z==1)
time1<-time[to.get]
delta1<-delta[to.get]

if(L.type=="quantile"){
L.0<-quantile(time0[delta0==1],0.975)
L.1<-quantile(time1[delta1==1],0.975)
L<-min(L.0,L.1)
}


if(L.type=="events"){
  L.0<-max(time0[delta0==1])
  L.1<-max(time1[delta1==1])
  L<-min(L.0,L.1)
}

if(L.type=="observed"){
  L.0<-max(time0)
  L.1<-max(time1)
  L<-min(L.0,L.1)
}

}


########################################################################################
#                           Re-define truncated versions                               #
# For estimation of truncated means (cf, Zhao & Tsiatis (2000) and Leon & Tsai (2004)) #
########################################################################################
delta.L<-1-(1-delta)*c(time<L)
time.L<-pmin(time,L)

################
# Control group
################
to.get<-which(z==0)
time0<-time.L[to.get]
delta0<-delta.L[to.get]
#####################
# Experimental group
####################
to.get<-which(z==1)
time1<-time.L[to.get]
delta1<-delta.L[to.get]

# Control KM & mean(L)
temp<-NA.CHR.Weighted.RMST(time=time0,delta=(delta0==1))
S0.KM<-temp$S.KM
se0.KM<-temp$se.KM
t0.points<-temp$at.points
R0<-temp$n.risk
dN0<-temp$n.events
# May need modification (R0-1 term) for ties
dC0<-ifelse(R0>1,dN0/(R0*(R0-1)),0.0)
###################################
# Restricted mean (point estimate)
###################################
m0.L<-LS.int(S0.KM,t0.points,L)

m0.t<-unlist(lapply(as.list(t0.points),mu.L,S=S0.KM,t=t0.points))
# Gill notation
mbar.0<-c(m0.L-m0.t)
var.m0L<-LS.int.del(fx=mbar.0^2,dx=dC0,x=t0.points,L=L)


# Experimental
temp<-NA.CHR.Weighted.RMST(time=time1,delta=(delta1==1))
S1.KM<-temp$S.KM
se1.KM<-temp$se.KM
t1.points<-temp$at.points
R1<-temp$n.risk
dN1<-temp$n.events
# May need modification (time1-1 term)for ties
dC1<-ifelse(R1>1,dN1/(R1*(R1-1)),0.0)
m1.L<-LS.int(S1.KM,t1.points,L)

m1.t<-unlist(lapply(as.list(t1.points),mu.L,S=S1.KM,t=t1.points))
mbar.1<-c(m1.L-m1.t)
var.m1L<-LS.int.del(fx=mbar.1^2,dx=dC1,x=t1.points,L=L)

dhat<-m1.L-m0.L

se.dhat.asy<-sqrt(var.m0L+var.m1L)
lower.asy<-dhat-1.96*se.dhat.asy
upper.asy<-dhat+1.96*se.dhat.asy

Z.asy<-dhat/se.dhat.asy
pval.asy.2sided<-2*(1-pnorm(abs(Z.asy)))
pval.asy<-1-pnorm(Z.asy)

if(details){
cat("Asymptotic based CI's","\n")
cat("Restricted means (dhat) at L=",c(L),"\n")
cat("dhat, Asy 95% CI=",c(dhat,lower.asy,upper.asy),"\n")
cat("Asy: Z, p-val (1-sided, 2-sided)=",c(Z.asy,pval.asy,pval.asy.2sided),"\n")
}

se.dhat<-se.dhat.asy

if(draws>0 | !is.null(G0.draws) | !is.null(G1.draws)){
# Resampling
# Control
temp<-KM.resample(time=time0,delta=(delta0==1),draws=draws,G.draws=G0.draws,details=details,at.points=t0.points)
#[1] "S.KM"      "at.points" "S.KM.star"
#W.1.star<-apply(S2.z0.star,2,surv.tau,tpoints=at.points,tau=0.5)
#se1.1<-sqrt(var(W.1.star))
# Distribution of (m0.hat-m0)
m.star.0<-apply(temp$S.KM.star,2,LS.int,x=temp$at.points,L=L)
# Experimental
temp<-KM.resample(time=time1,delta=(delta1==1),draws=draws,G.draws=G1.draws,details=details,at.points=t1.points)
# Distribution of (m1.hat-m1)
m.star.1<-apply(temp$S.KM.star,2,LS.int,x=temp$at.points,L=L)
# 95% CIs
dhat.star<-m.star.1-m.star.0
Zdhat.star<-dhat.star/se.dhat.asy
se.dhat.star<-sqrt(var(dhat.star))
# Normal-type
lower.star<-dhat-1.96*sqrt(var(dhat.star))
upper.star<-dhat+1.96*sqrt(var(dhat.star))
if(details){ 
cat("Martingale-resampling based CI's","\n")  
cat("Number of resampling draws=",c(draws),"\n")
cat("Restricted means (dhat) at L=",c(L),"\n")
cat("dhat, 95% CI=",c(dhat,lower.star,upper.star),"\n")
}
# Note: Can also calculate percentile versions of CIs
}

if(get.band){
# Non-truncated version

################
# Control group
################
to.get<-which(z==0)
time0<-time[to.get]
delta0<-delta[to.get]
#####################
# Experimental group
####################
to.get<-which(z==1)
time1<-time[to.get]
delta1<-delta[to.get]

if(is.null(taus)){
tau1.0<-quantile(time0[delta0==1],0.025)
tau1.1<-quantile(time1[delta1==1],0.025)
tau1<-max(c(tau1.0,tau1.1))

tau2.0<-quantile(time0[delta0==1],0.975)
tau2.1<-quantile(time1[delta1==1],0.975)
tau2<-min(c(tau2.0,tau2.1))
}

if(!is.null(taus)){
tau1<-taus[1]
tau2<-taus[2]
}

if(details) cat("Resampling for simultaneous CIs over (tau1,tau2)","\n")

# Points at which to calculated difference in KMs
dpoints<-seq(tau1,tau2,by=tau.seq)

temp<-NA.CHR.Weighted.RMST(time=time0,delta=(delta0==1),at.points=dpoints)
S0.dpoints<-temp$S.KM
se0.dpoints<-temp$se.KM

temp<-NA.CHR.Weighted.RMST(time=time1,delta=(delta1==1),at.points=dpoints)
S1.dpoints<-temp$S.KM
se1.dpoints<-temp$se.KM

diff.dpoints<-S1.dpoints-S0.dpoints
sig2.dpoints<-se1.dpoints^2+se0.dpoints^2

temp<-KM.resample(time=time0,delta=(delta0==1),draws=draws.band,at.points=dpoints,details=details)
S0.star.dpoints<-temp$S.KM.star

temp<-KM.resample(time=time1,delta=(delta1==1),draws=draws.band,at.points=dpoints,details=details)
S1.star.dpoints<-temp$S.KM.star

Z.star.dpoints<-(S1.star.dpoints-S0.star.dpoints)/sqrt(sig2.dpoints)

sups<-apply(abs(Z.star.dpoints),2,max,na.rm=TRUE)
calpha<-quantile(sups,c(0.95))

# pointwise CIs
pw.lower<-diff.dpoints-1.96*sqrt(sig2.dpoints)
pw.upper<-diff.dpoints+1.96*sqrt(sig2.dpoints)

# simulataneous band
sb.lower<-diff.dpoints-calpha*sqrt(sig2.dpoints)
sb.upper<-diff.dpoints+calpha*sqrt(sig2.dpoints)

if(plotband){
plot.band(x=dpoints,mean=diff.dpoints,xlabel="time(t)",
ylabel=expression(delta(t)==hat(S)[1](t)-hat(S)[0](t)),
band=TRUE,lower=sb.lower,upper=sb.upper,show.axes=TRUE,ltype="l")
lines(dpoints,pw.upper,type="l",lty=2,lwd=0.50)
lines(dpoints,pw.lower,type="l",lty=2,lwd=0.50)
abline(h=0.0,lty=1,col="black",lwd=2)
rug(dpoints)
}

if(details) cat("Simultaneous confidence band over [tau1,tau2]=",c(tau1,tau2),"\n")

}

if(draws>0 & get.band==FALSE) return(list(dhat=dhat,m1L=m1.L,m0L=m0.L,m0.star=m.star.0,m1.star=m.star.1,dstar=dhat.star,Zdstar=Zdhat.star,
                                          lower.asy=lower.asy,upper.asy=upper.asy,
                                          lower.star=lower.star, upper.star=upper.star,
                                          L=L,Z.asy=Z.asy,pval.asy=pval.asy,se.dhat.star=se.dhat.star,se.dhat=se.dhat))
if(draws==0 & get.band==FALSE) return(list(dhat=dhat,m1L=m1.L,m0L=m0.L,L=L,Z.asy=Z.asy,pval.asy=pval.asy,se.dhat=se.dhat,lower=lower.asy,upper=upper.asy))

if(draws>0 & get.band==TRUE) return(list(dhat=dhat,m1L=m1.L,m0L=m0.L,m0.star=m.star.0,m1.star=m.star.1,dstar=dhat.star,Zdstar=Zdhat.star,
lower.asy=lower.asy,upper.asy=upper.asy,lower.star=lower.star, upper.star=upper.star,L=L,
dpoints=dpoints,S1.dpoints=S1.dpoints,S0.dpoints=S0.dpoints,pw.lower=pw.lower,pw.upper=pw.upper,sb.lower=sb.lower,sb.upper=sb.upper,tau1=tau1,tau2=tau2,calpha=calpha,Z.star.dpoints=Z.star.dpoints,sup.star=sups,
diff.dpoints=diff.dpoints,sig2.dpoints=sig2.dpoints,Z.asy=Z.asy,pval.asy=pval.asy,se.dhat.star=se.dhat.star,se.dhat=se.dhat))

if(draws==0 & get.band==TRUE) return(list(dhat=dhat,m1L=m1.L,m0L=m0.L,L=L,lower=lower.asy,upper=upper.asy,
dpoints=dpoints,S1.dpoints=S1.dpoints,S0.dpoints=S0.dpoints,pw.lower=pw.lower,pw.upper=pw.upper,sb.lower=sb.lower,sb.upper=sb.upper,tau1=tau1,tau2=tau2,calpha=calpha,Z.star.dpoints=Z.star.dpoints,sup.star=sups,
diff.dpoints=diff.dpoints,sig2.dpoints=sig2.dpoints,Z.asy=Z.asy,pval.asy=pval.asy,se.dhat=se.dhat))
}
