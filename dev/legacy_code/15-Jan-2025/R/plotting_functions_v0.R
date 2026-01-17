
N.Weighted<-function(x,error,W=rep(1,length(error))){
  sum(W*(error<=x))
}

R.Weighted<-function(x,error,W=rep(1,length(error))){
  sum(W*(error>=x))
}

NA.CHR.Weighted<-function(time,Delta,W.n=rep(1,length(time)),W.d=rep(1,length(time)),
                          at.points=sort(time),se.type="greenwood",get.Stute=FALSE,tpoints.add=NULL){
  
  if(!is.null(tpoints.add)) at.points<-sort(c(unique(c(at.points,tpoints.add))))
  
  if(se.type!="greenwood" & se.type!="tsiatis") stop("Invalid se type -- greenwood or tsiatis allowed")
  #is.sorted<-(all(time==sort(time)))
  is.sorted<-!is.unsorted(time)
  if(!is.sorted){
    id<-order(time); time<-time[id]; Delta<-Delta[id]; W.n<-W.n[id]; W.d<-W.d[id]
  }
  risk<-unlist(lapply(as.list(at.points),R.Weighted,error=time,W=W.d))
  ###########################################################################
  ### Adaptive H process correspoinding to N-A rep. via integral wrt M.G ####
  Hmart.chf<-ifelse(risk>0,1/risk,0)
  ############################################################################
  counting<-unlist(lapply(as.list(at.points),N.Weighted,error=time,W=W.n*ifelse(Delta==1,1,0)))
  counting <- c(0, counting)
  dN<-diff(counting)
  dN.risk<-ifelse(risk>0,dN/risk,0.0)
  chf <- cumsum(dN.risk)
  var.chf<-cumsum(ifelse(risk>0,dN/(risk^2),0.0))
  S.KM <- cumprod(1-dN.risk)
  
  S.KM[which(S.KM<0)]<-0.0
  
  S.NA <- exp(-chf)
  var.NA<-(S.NA^2)*var.chf
  # Greenwood variance estimate
  if(se.type=="greenwood"){
    aa<-dN
    bb<-risk*(risk-dN)
    var.KM<-(S.KM^2)*cumsum(ifelse(risk>0,aa/bb,0.0))
    se.KM<-sqrt(var.KM)
  }
  if(se.type=="tsiatis"){
    var.KM<-(S.KM^2)*var.chf
    se.KM<-sqrt(var.KM)
  }
  result<-list(time=time,at.points=at.points,S.NA=S.NA,S.KM=S.KM,chf=chf,se.chf=sqrt(var.chf),
               se.NA=sqrt(var.NA),dN.risk=dN.risk,
               n.risk=risk,dN=dN,
               Hmart.chf=Hmart.chf,se.KM=se.KM)
  return(result)
}

# Some plotting functions

plot_twosample<-function(df,tte.name,event.name,treat.name,wgt.name=NULL,xloc1=NULL,xloc2=NULL,details=FALSE,show.logrank=FALSE,
                         col.treat="blue",col.control="grey",choose_ylim=FALSE,
                         show.cox=TRUE,censor.cex=1,cox.cex=0.7,prob.points=c(0,0.25,0.5,0.75,1.0),
                         show.Y.axis=TRUE, cex_Yaxis=1,
                         exp.lab="Treat",con.lab="Control",legend.cex=0.70,risk.cex=0.65,yloc1=0.6,yloc2=0.6,subid=NULL,byrisk=2,fix.rows=TRUE,
                         show.med=TRUE,med.cex=1,med.offset=0.25,
                         ylab="Survival",xlab="Months"){
  if(is.null(wgt.name)){
    df$wgt<-rep(1,nrow(df)) 
  }
  tpoints.add<-c(-1)
  km.fit<-KM.plot.2sample.weighted(Y=df[,c(tte.name)],E=df[,c(event.name)],Treat=df[,c(treat.name)],Weight=df$wgt,
                                   col.2=col.control,col.1=col.treat,choose_ylim=choose_ylim,
                                   risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
                                   prob.points=prob.points,cex_Yaxis=cex_Yaxis,show.Y.axis=show.Y.axis,
                                   stop.onerror=TRUE,Xlab=xlab,Ylab=ylab,details=details,censor.cex=censor.cex,
                                   show.logrank=show.logrank,show.med=FALSE,med.cex=med.cex,med.offset=med.offset,
                                   show.cox=show.cox,cox.cex=cox.cex)
  cpoints <- km.fit$cpoints
  
  m1<-round(km.fit$med.1,2)
  m0<-round(km.fit$med.2,2)
  
  if(is.null(xloc1)) xloc1<-c(diff(range(cpoints))/2)  
  
  if(show.med){
  legend("top", c(paste(c(exp.lab,eval(parse(text=m1))),collapse="; med="), paste(c(con.lab,eval(parse(text=m0))),collapse="; med=")), 
  col = c(col.treat,col.control), lwd = 2, bty = "n", cex=legend.cex)
  }
  if(!show.med) legend("top", c(exp.lab, con.lab), col = c(col.treat,col.control), lwd = 2, bty = "n", cex=legend.cex)
  
}  


# For plotting subgroups and the complement

plot.subgroup<-function(tte.name,event.name,treat.name,wgt.name=NULL,sub1,sub1C,xloc1=NULL,xloc2=NULL,details=FALSE,show.logrank=FALSE,
                        ymin=0,cox.cex=0.6,prob.points=c(0,0.25,0.5,0.75,1.0),cex_Yaxis=1,title1=NULL,title2=NULL,choose_ylim=TRUE,
                        exp.lab="Treat",con.lab="Control",legend.cex=0.60,risk.cex=0.60,yloc1=0.6,yloc2=0.6,subid=NULL,byrisk=2,fix.rows=TRUE,show.med=FALSE,ylab="Survival"){
  
  
  if(!is.null(subid) & (!is.null(title1) | !is.null(title2))) stop("subid OR titles need to be specified (not both)")
  
  if(is.null(wgt.name)){
    sub1$wgt<-rep(1,nrow(sub1)) 
    sub1C$wgt<-rep(1,nrow(sub1C))
  }
  
  if(!is.null(wgt.name)){
    sub1$wgt<-sub1[,c(wgt.name)] 
    sub1C$wgt<-sub1C[,c(wgt.name)]
  }
  
  # Set tpoints.add to span range of both subgroups
  aa<-max(sub1[,c(tte.name)])
  bb<-max(sub1C[,c(tte.name)])
  tpoints.add<-c(-1,max(aa,bb))
  
  if(fix.rows) par(mfrow=c(1,2))
  df<-sub1
 
   km.fit<-KM.plot.2sample.weighted(Y=df[,c(tte.name)],E=df[,c(event.name)],Treat=df[,c(treat.name)],Weight=df$wgt,
                                   risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
                                   cex_Yaxis=cex_Yaxis,choose_ylim=choose_ylim,
                                   stop.onerror=TRUE,Xlab="Months",Ylab=ylab,details=details,prob.points=prob.points,
                                   ymin=ymin,show.logrank=show.logrank,show.med=TRUE,show.cox=TRUE,cox.cex=cox.cex)
  cpoints <- km.fit$cpoints
  
  m1<-round(km.fit$med.1,1)
  m0<-round(km.fit$med.2,1)
  
  if(is.null(xloc1)) xloc1<-c(diff(range(cpoints))/2)  
  
  #if(!show.med) legend(xloc1, yloc1, c(paste(c(exp.lab,eval(parse(text=m1))),collapse="; med="), paste(c(con.lab,eval(parse(text=m0))),collapse="; med=")), 
  #                     col = c("black","blue"), lwd = 2, bty = "n", cex=legend.cex)
  #if(show.med) legend(xloc1, yloc1, c(exp.lab, con.lab), col = c("black","blue"), lwd = 2, bty = "n")
  
  if(show.med){
  legend("top", c(paste(c(exp.lab,eval(parse(text=m1))),collapse="; med="), paste(c(con.lab,eval(parse(text=m0))),collapse="; med=")), 
  col = c("black","blue"), lwd = 2, bty = "n", cex=legend.cex)
  }
  if(!show.med) legend("top", c(exp.lab, con.lab), col = c("black","blue"), lwd = 2, bty = "n", cex=legend.cex)
  
  #legend("top",c(subid),bty="n",cex=0.7)
  if(!is.null(subid)) title(main=subid,cex.main=0.80)
  if(!is.null(title1)) title(main=title1,cex.main=0.80)
  
  df<-sub1C
  
  km.fit<-KM.plot.2sample.weighted(Y=df[,c(tte.name)],E=df[,c(event.name)],Treat=df[,c(treat.name)],Weight=df$wgt,
                                   risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
                                   stop.onerror=TRUE,Xlab="Months",Ylab="",details=details,show.Y.axis=FALSE,
                                   ymin=ymin,cox.cex=cox.cex,cex_Yaxis=cex_Yaxis,choose_ylim=choose_ylim,
                                   show.logrank=show.logrank,show.med=TRUE,show.cox=TRUE)
  
  cpoints <- km.fit$cpoints
  m1<-round(km.fit$med.1,1)
  m0<-round(km.fit$med.2,1)
  
  if(is.null(xloc2)) xloc2<-c(diff(range(cpoints))/2)  
  
  if(show.med){
  legend("top", c(paste(c(exp.lab,eval(parse(text=m1))),collapse="; med="), paste(c(con.lab,eval(parse(text=m0))),collapse="; med=")), 
   col = c("black","blue"), lwd = 2, bty = "n",cex=legend.cex)
  }  
  if(!show.med) legend("top", c(exp.lab, con.lab), col = c("black","blue"), lwd = 2, bty = "n", cex=legend.cex)
  
  
  #if(!show.med) legend(xloc2, yloc2, c(paste(c(exp.lab,eval(parse(text=m1))),collapse="; med="), paste(c(con.lab,eval(parse(text=m0))),collapse="; med=")), 
  #                     col = c("black","blue"), lwd = 2, bty = "n",cex=legend.cex)
  #if(show.med) legend(xloc2, yloc2, c(exp.lab, con.lab), col = c("black","blue"), lwd = 2, bty = "n")
  
  
  if(!is.null(subid)) title(main="Complement",cex.main=0.80)
  if(!is.null(title2)) title(main=title2,cex.main=0.80)
}

# KM plot functions
###################################################################################################################################################
# K-M is weighted via Weight (Weight==1 is standard KM)
# Cox model is either via Weighting OR Covariates (if Covariates = NULL, then Cox model is weighted;  Otherwise, Cox model included Covariates)
###################################################################################################################################################

KM.plot.2sample.weighted<-function(Y,E,Treat,Covariates=NULL,Weight=rep(1,length(Y)),show.cox=FALSE,col.cox="red",cox.cex=0.725,show.logrank=FALSE,stop.onerror=FALSE,check.KM=FALSE,
                                   details=FALSE,put.legend.cox='topright',put.legend.lr='top',del.med=0.075,lr.digits=2,cox.digits=3,
                                   tpoints.add=c(-1),by.risk=NULL,Xlab="time",Ylab="proportion surviving",col.1="black",col.2="blue",show.med=FALSE,
                                   choose_ylim =FALSE,
                                   med.cex=1,med.offset=0.25,risk.cex=1,
                                   plotit=TRUE,
                                   ltys=c(1,1),lwds=c(1,1),
                                   quant=0.5,
                                   censor.mark.all=TRUE,censor.cex=1.0,
                                   show.ticks=TRUE,risk.set=TRUE,
                                   ymin=0,ymax=1,
                                   ymin2=NULL, risk_offset=0.075, risk_delta=0.025,
                                   y.risk2=NULL,show.Y.axis=TRUE,cex_Yaxis=1,
                                   y.risk1=NULL,
                                   add.segment=FALSE,risk.add=NULL,xmin=0,xmax=NULL,x.truncate=NULL,time.zero=0.0,prob.points=c(seq(0.1,1.0,by=0.1))){
  
  #if(censor.mark.all==TRUE) warning("All censorings will be marked regardless of *tying to an event* [DOES NOT REPLICATE R Survfit()]")
  Weight_orig <- Weight
  Y1<-Y[which(Treat==1)]
  E1<-E[which(Treat==1)]
  W1<-Weight[which(Treat==1)]
  
  Y0<-Y[which(Treat==0)]
  E0<-E[which(Treat==0)]
  W0<-Weight[which(Treat==0)]
  
  Time<-c(Y1,Y0)
  Event<-c(E1,E0)
  Strata<-c(rep(2,length(Y1)),rep(1,length(Y0)))
  Weight<-c(W1,W0)
  if(!is.null(Covariates)){ 
    coxfit<-coxph(Surv(Y,E)~Treat+Covariates)  
    hr.ci<-as.matrix(round(summary(coxfit)$conf.int,3))
    pval.cox<-summary(coxfit)$coefficients[1,6]
  }
  if(is.null(Covariates)){ 
    coxfit<-coxph(Surv(Y,E)~Treat)  
    hr.ci<-as.matrix(round(summary(coxfit)$conf.int,3))
    pval.cox<-summary(coxfit)$coefficients[5]
  }
  # ipw weighting
  if(is.null(Covariates) & !all(Weight_orig==1)){ 
  coxfit <- coxph(Surv(Y,E)~Treat, weights=Weight_orig, robust=TRUE)
  hr.ci<-as.matrix(round(summary(coxfit)$conf.int,3))
  pval.cox<-summary(coxfit)$coefficients[6]
    }
  
  # Get HR for Treat
  hr.ci<-hr.ci[1,c(1,3,4)]
  
  if(details){
    if(!is.null(Covariates)){
      cat("Cox adjusted HR=",c(paste(hr.ci[1])),"\n")
      cat("Cox CIs=",c(paste(hr.ci[2]),paste(hr.ci[3])),"\n")}
      if(is.null(Covariates) & all(Weight_orig==1)){
      cat("Cox un-adjusted HR=",c(paste(hr.ci[1])),"\n")
      cat("Cox CIs=",c(paste(hr.ci[2]),paste(hr.ci[3])),"\n")}
    
    if(is.null(Covariates) & !all(Weight_orig==1)){
      cat("Cox ipw-adjusted HR=",c(paste(hr.ci[1])),"\n")
      cat("Cox CIs=",c(paste(hr.ci[2]),paste(hr.ci[3])),"\n")}
  }
  
  
  if(is.null(by.risk)){
    tt<-sort(unique(Time))
    by.risk<-quantile(tt-c(0,tt[-length(tt)]),c(0.75),na.rm=TRUE)
  }
  
  stratums<-c(sort(unique(Strata),decreasing=TRUE))
  if(any(is.na(Strata))) stop("Missing strata information")
  
  x1<-Time[which(Strata==stratums[1])]
  to.get<-!is.na(x1)
  x1<-x1[to.get]
  eta1<-Event[which(Strata==stratums[1])]
  eta1<-eta1[to.get]
  w1<-Weight[which(Strata==stratums[1])]
  w1<-w1[to.get]
  
  x2<-Time[which(Strata==stratums[2])]
  to.get<-!is.na(x2)
  x2<-x2[to.get]
  eta2<-Event[which(Strata==stratums[2])]
  eta2<-eta2[to.get]
  w2<-Weight[which(Strata==stratums[2])]
  w2<-w2[to.get]
  
  # Common points to estimate both arms to include all event times (for differences)
  tpoints.add<-sort(unique(c(tpoints.add,sort(unique(Y[which(E==1)])))))
  
  xx<-c(x1,x2)
  ee<-c(eta1,eta2)
  at.points<-sort(c(unique(c(x1,x2,tpoints.add))))
  at.points1<-sort(c(unique(c(x1,tpoints.add))))
  at.points2<-sort(c(unique(c(x2,tpoints.add))))
    # Experimental
  fit1<-NA.CHR.Weighted(time=x1,Delta=(eta1==1),W.n=w1,W.d=w1,at.points=at.points1)
  S1.KM<-fit1$S.KM
  SE1.KM<-fit1$se.KM
  
  # Control
  fit2<-NA.CHR.Weighted(time=x2,Delta=(eta2==1),W.n=w2,W.d=w2,at.points=at.points2)
  S2.KM<-fit2$S.KM
  SE2.KM<-fit2$se.KM
  
  if(choose_ylim){
  ymin <- min(c(S1.KM,S2.KM))
  }
  
  if(is.null(ymin2)) ymin2 <- ymin-risk_offset  
  if(is.null(y.risk1)) y.risk1 <- ymin2
  if(is.null(y.risk2)) y.risk2 <- ymin2+risk_delta
  
  # Points for plotting
  #tpoints.add.plot <- tpoints.add
  # Same at at.points isf default
  at.points1.plot<-sort(c(unique(c(x1,tpoints.add))))
  at.points2.plot<-sort(c(unique(c(x2,tpoints.add))))
  
  # LR() denotes the log-rank scale
  # Integral of "LR differences" = log-rank test (non-standardized)
  ZLR.dpoints<-LR.dpoints<-S1.dpoints<-SE1.dpoints<-S2.dpoints<-SE2.dpoints<-NULL
  
  tp1.match<-match(tpoints.add,at.points1)
  tp2.match<-match(tpoints.add,at.points2)
  
  S1.dpoints<-S1.KM[tp1.match]
  SE1.dpoints<-SE1.KM[tp1.match]
  Y1.dpoints<-fit1$n.risk[tp1.match]
  dN1.dpoints<-fit1$dN[tp1.match]
  
  S2.dpoints<-S2.KM[tp2.match]
  SE2.dpoints<-SE2.KM[tp2.match]
  Y2.dpoints<-fit2$n.risk[tp2.match]
  dN2.dpoints<-fit2$dN[tp2.match]
  
  K<-(Y1.dpoints*Y2.dpoints)/(Y1.dpoints+Y2.dpoints)
  
  Y.dpoints<-Y1.dpoints+Y2.dpoints
  dN.dpoints<-dN1.dpoints+dN2.dpoints
  
  term1<-cumsum(ifelse(Y1.dpoints>0,(K/Y1.dpoints)*dN1.dpoints,0.0))
  term2<-cumsum(ifelse(Y2.dpoints>0,(K/Y2.dpoints)*dN2.dpoints,0.0))
  
  # variance
  h0<-ifelse(Y1.dpoints==0,0,(K^2/Y1.dpoints))
  h1<-ifelse(Y2.dpoints==0,0,(K^2/Y2.dpoints))
  dJ<-ifelse(Y.dpoints==1,0,(dN.dpoints-1)/(Y.dpoints-1))
  dL<-ifelse(Y.dpoints==0,0,dN.dpoints/Y.dpoints)
  sig2s<-(h0+h1)*(1-dJ)*dL
  sig2<-cumsum(sig2s)
  
  LR.dpoints<-term2-term1
  ZLR.dpoints<-LR.dpoints/sqrt(sig2)
  # Zlr() is weighted log-rank process
  # Zlr(tau) for tau end of followup is 
  # the "final" z statistic (1-sided)
  z.lr<-ZLR.dpoints[length(ZLR.dpoints)]
  p.lr<-1-pnorm(z.lr)
  
  if(check.KM & all(Weight==1)){
    p.2side <- 1-pchisq(z.lr^2,1)
    logrnk<-survdiff(Surv(Y,E)~Treat)
    p.val<-1-pchisq(logrnk$chisq,length(logrnk$n) - 1)
    if(round(p.2side-p.val,8)>0){
    warning("Discrepancy with log-rank process and survdiff")
    if(details){
    cat("My p-value and survdiff=",c(p.2side,p.val),"\n")  
    cat("My z^2 and survdiff=",c(z.lr^2,logrnk$chisq),"\n")
    }    
    }  
  }
  
  
  #########################################
  # Note: want to confirm with R function
  #########################################
  #           Check KM fits               #
  if(check.KM){
    # only at observed unique observations to align with survfit
    # Experimental
    S1.check<-NA.CHR.Weighted(time=x1,Delta=(eta1==1),W.n=w1,W.d=w1,at.points=sort(unique(x1)))$S.KM
    S2.check<-NA.CHR.Weighted(time=x2,Delta=(eta2==1),W.n=w2,W.d=w2,at.points=sort(unique(x2)))$S.KM
    KM.fit<-survfit(Surv(Y1,E1)~1,weights=W1)
    max.error<-round(max(abs(KM.fit$surv-S1.check)),8)
    if(details) print(max.error)
    if(stop.onerror & max.error>=0.00001) stop("Discrepancy in KM1 fit")
    KM.fit<-survfit(Surv(Y0,E0)~1,weights=W0)
    max.error<-round(max(abs(KM.fit$surv-S2.check)),8)
    if(details) print(max.error)
    if(stop.onerror & max.error>=0.00001) stop("Discrepancy in KM2 fit")
    #        End Check KM fits               #
  }
  events1<-sort(unique(x1[eta1==1]))
  events2<-sort(unique(x2[eta2==1]))
  
  if(!censor.mark.all){
    # Censored points
    # Curves are marked at each censoring time which is not also a death time
    x1.cens<-x1[eta1==0] 
    x2.cens<-x2[eta2==0]
    x1.cens<-setdiff(x1.cens,events1)  # Censorings that are NOT an event (this is the R version which follows estimation convention)
    x2.cens<-setdiff(x2.cens,events2)
    x1.match<-match(x1.cens,at.points1)
    x2.match<-match(x2.cens,at.points2)
  }
  
  if(censor.mark.all){
    # Censored points
    x1.cens<-x1[eta1==0] 
    x2.cens<-x2[eta2==0]
    x1.match<-match(x1.cens,at.points1)
    x2.match<-match(x2.cens,at.points2)
  }
  
  risk.points<-round(c(seq(0,max(at.points),by=by.risk)))
  risk.points<-sort(unique(c(risk.points,risk.add)))
  
  risk.2<-unlist(lapply(as.list(risk.points),R.Weighted,error=x2,W=w2))
  risk.1<-unlist(lapply(as.list(risk.points),R.Weighted,error=x1,W=w1))
  
  risk.2<-round(risk.2)
  risk.1<-round(risk.1)
  
  # Note: "med" denotes general quantile, but median is default
  med.1<-suppressWarnings(min(at.points1[S1.KM<=quant]))
  # log transform
  KM1.lower<-exp(log(S1.KM)-1.96*SE1.KM/S1.KM)
  KM1.upper<-exp(log(S1.KM)+1.96*SE1.KM/S1.KM)
  med.1.lower<-suppressWarnings(min(at.points1[KM1.lower<=quant]))
  med.1.upper<-suppressWarnings(min(at.points1[KM1.upper<=quant]))
  
  if(details){
    cat("Stratum,n,events=",c(stratums[1],length(x1),sum(eta1)),"\n")
    cat("Median, Lower, Upper=",c(med.1,med.1.lower,med.1.upper),"\n")
  }
  
  med.2<-suppressWarnings(min(at.points2[S2.KM<=quant]))
  KM2.lower<-exp(log(S2.KM)-1.96*SE2.KM/S2.KM)
  KM2.upper<-exp(log(S2.KM)+1.96*SE2.KM/S2.KM)
  med.2.lower<-suppressWarnings(min(at.points2[KM2.lower<=quant]))
  med.2.upper<-suppressWarnings(min(at.points2[KM2.upper<=quant]))
  
  if(details){
    cat("Stratum,n,events=",c(stratums[2],length(x2),sum(eta2)),"\n")
    cat("Median, Lower, Upper=",c(med.2,med.2.lower,med.2.upper),"\n")
  }
  # Experimental
  fit1<-NA.CHR.Weighted(time=x1,Delta=(eta1==1),W.n=w1,W.d=w1,at.points=at.points1.plot)
  S1.KM.plot<-fit1$S.KM
  # Control
  fit2<-NA.CHR.Weighted(time=x2,Delta=(eta2==1),W.n=w2,W.d=w2,at.points=at.points2.plot)
  S2.KM.plot<-fit2$S.KM
  Y1.plot<-S1.KM.plot 
 Y2.plot<-S2.KM.plot
  
  Y1<-S1.KM
  Y2<-S2.KM
  
  if(is.null(x.truncate)) xmax<-max(c(at.points,xmax))
  if(!is.null(x.truncate)) xmax<-x.truncate
  
  if(plotit){
    plot(at.points1.plot,Y1.plot,type="s",ylim=c(ymin2,ymax),xlim=c(xmin,xmax),lty=ltys[1],col=col.1,lwd=lwds[1],xlab=Xlab,ylab=Ylab,axes=FALSE)
    lines(at.points2.plot,Y2.plot,lty=ltys[2],type="s",col=col.2,lwd=lwds[2])
    
    if(show.ticks==TRUE){
      points(x1.cens,Y1[x1.match],pch=3,col=col.1,cex=censor.cex)
      points(x2.cens,Y2[x2.match],pch=3,col=col.2,cex=censor.cex)
    }
    
    #if(what.toplot=="KM") abline(v=0,col="grey",lwd=3,lty=2)
    
    abline(h=ymin-0.035,lty=1,col=1)
    #axis(2,at=c(0.0,0.2,0.4,0.5,0.6,0.8,1.0))
    
    # Include at least 2-points in-between ymin and ymax
    if((ymax-ymin2) <= 0.5) prob.points <- round(seq(ymin,ymax,length=4),2)
    
    if(show.Y.axis) axis(2,at=c(prob.points),cex.axis=cex_Yaxis)
    axis(1,at=risk.points,labels=c(risk.points+time.zero),cex.axis=cex_Yaxis)
    #axTicks(side=1, axp =c(risk.points+5.25))
    
    box()
    if(risk.set){
      text(c(risk.points),c(y.risk2),c(risk.2),col=col.2,cex=risk.cex)
      text(c(risk.points),c(y.risk1),c(risk.1),col=col.1,cex=risk.cex)
    }
    if(add.segment) segments(-10,1,0,1)  
    
    if(show.med){
      if(med.1!=Inf & med.2!=Inf){
        text(xmax-4,ymax-med.offset,paste(round(med.1,1)),col=col.1,cex=med.cex)
        text(xmax-4,ymax-(med.offset+del.med),paste(round(med.2,1)),col=col.2,cex=med.cex)
      }
#     if(med.1!=Inf & med.2==Inf){
#        text(xmax-2,ymax-0.2,paste(round(med.1,1)),col=col.1,cex=med.cex)
#        text(xmax-2,ymax-(0.2+del.med),paste("median NA"),col=col.2,cex=1)
#      }
#      if(med.1==Inf & med.2!=Inf){
#        text(xmax-2,ymax-0.2,paste("median NA"),col=col.1,cex=1)
#        text(xmax-2,ymax-(0.2+del.med),paste(round(med.2,1)),col=col.2,cex=med.cex)
#      }
#      if(med.1==Inf & med.2==Inf){
#        text(xmax-4,ymax-0.2,paste("median NA"),col=col.1,cex=1)
#        text(xmax-4,ymax-(0.2+del.med),paste("median NA"),col=col.2,cex=1)
#      }
      
    }
    if(show.cox){
      rp<-vector('expression',4)
      rp[1]=substitute(expression(italic(hr)== MYVALUE1), 
                       list(MYVALUE1 = format(hr.ci[1],digits = cox.digits)))[2]
      rp[2]=substitute(expression(italic(lb) == MYVALUE2), 
                       list(MYVALUE2 = format(hr.ci[2], digits = cox.digits)))[2]
      rp[3]=substitute(expression(italic(ub) == MYVALUE3), 
                       list(MYVALUE3= format(hr.ci[3], digits = cox.digits)))[2]               
      rp[4]=substitute(expression(italic(p) == MYVALUE3), 
                       list(MYVALUE3= format(pval.cox, digits = cox.digits)))[2]               
      legend(put.legend.cox, legend = rp, bty = 'n',cex=cox.cex)
    }
    
    if(show.logrank){
      rp<-vector('expression',1)
      rp[1]=substitute(expression(italic(p[superiority])== MYVALUE), 
                       list(MYVALUE = format(p.lr,digits=lr.digits)))[2]
      legend(put.legend.lr, legend = rp, bty = 'n',cex=cox.cex)
    }
  }
  
  return(list(S1.KM=S1.KM,S2.KM=S2.KM,KM1.lower=KM1.lower,KM2.lower=KM2.lower,KM1.upper=KM1.upper,KM2.upper=KM2.upper,
              pval.cox=pval.cox,pval.lr=p.lr,cpoints=tpoints.add,
              S1.dpoints=S1.dpoints,SE1.dpoints=SE1.dpoints,
              S2.dpoints=S2.dpoints,SE2.dpoints=SE2.dpoints,
              LR.dpoints=LR.dpoints,ZLR.dpoints=ZLR.dpoints,
              at.points1=at.points1,at.points2=at.points2,
              Y1.plot=Y1.plot,Y2.plot=Y2.plot,
              at.points1.plot=at.points1.plot,at.points2.plot=at.points2.plot,
              med.1=med.1,med.2=med.2,events.1=sum(eta1),events.2=sum(eta2),
              y.risk1=y.risk1,y.risk2=y.risk2,
              risk.1=risk.1,risk.2=risk.2,risk.points=risk.points))
}




hr.overlap<-function(data,type="histogram",stat.name="hr",strata.name="method", bins = 50L, 
                     alpha=0.5,colorit="grey50",titleit,xtitle="Hazard ratio estimates"){
  hr<-as.vector(data[,c(stat.name)])  
  method<-data[,c(strata.name)]
  hr.methods<-data.frame(Method = as.factor(method),hr.est = hr)
  Method<- hr.est <- NULL
  
  if (type == "density") {
    
    pl.obj <- ggplot(hr.methods, aes(x = hr.est, fill = Method)) + 
      geom_density(alpha = alpha, colour = colorit) + 
      geom_rug(aes(colour = Method)) + theme(legend.position = "bottom") + 
      ggtitle(titleit) + 
      xlab(xtitle)
    
  }
  else if (type == "histogram") {
    pl.obj <- ggplot(hr.methods, aes(x = hr.est, fill = Method)) + 
      geom_histogram(bins = bins, alpha = alpha, position = "identity") + 
      geom_rug(aes(colour = Method)) + theme(legend.position = "bottom") + 
      ggtitle("Histograms of estimated hrs") + 
      xlab("Hazard ratio estimates")
  }
  
  pl.obj    
}

or.overlap<-function(data,type="histogram",stat.name="or",strata.name="method", bins = 50L, alpha=0.5,colorit="grey50",titleit,vref=NULL){
  hr<-as.vector(data[,c(stat.name)])  
  method<-data[,c(strata.name)]
  hr.methods<-data.frame(Method = as.factor(method),hr.est = hr)
  Method<- hr.est <- NULL
  
  if (type == "density") {
    
    pl.obj <- ggplot(hr.methods, aes(x = hr.est, fill = Method)) + 
      geom_density(alpha = alpha, colour = colorit) + 
      geom_rug(aes(colour = Method)) + theme(legend.position = "bottom") + 
      ggtitle(titleit) + 
      xlab("Odds ratio estimates")
    
  }
  else if (type == "histogram") {
    pl.obj <- ggplot(hr.methods, aes(x = hr.est, fill = Method)) + 
      geom_histogram(bins = bins, alpha = alpha, position = "identity") + 
      geom_vline(xintercept = vref) +
      geom_rug(aes(colour = Method)) + theme(legend.position = "bottom") + 
      ggtitle("Histograms of estimated Odds ratios") + 
      xlab("Odds ratio estimates")
  }
  
  pl.obj    
}


hr.overlap2<-function(data,type="histogram",stat.name="hr",strata.name="method", bins = 50L, alpha=0.5,colorit="grey50",
                      titleit,xmin=0,xmax=20,vref=NULL){
  hr<-as.vector(data[,c(stat.name)])  
  method<-data[,c(strata.name)]
  hr.methods<-data.frame(Method = as.factor(method),hr.est = hr)
  Method<- hr.est <- NULL
  
  if (type == "density") {
    
    pl.obj <- ggplot(hr.methods, aes(x = hr.est, fill = Method)) + 
      geom_density(alpha = alpha, colour = colorit) + 
      xlim(xmin,xmax) +
      geom_vline(xintercept = vref) +
      geom_rug(aes(colour = Method)) + theme(legend.position = "bottom") + 
      ggtitle(titleit) + 
      xlab("Hazard ratio estimates")
    
  }
  else if (type == "histogram") {
    pl.obj <- ggplot(hr.methods, aes(x = hr.est, fill = Method)) + 
      geom_histogram(bins = bins, alpha = alpha, position = "identity") + 
      xlim(xmin,xmax) +
      geom_vline(xintercept = vref) +
      geom_rug(aes(colour = Method)) + theme(legend.position = "bottom") + 
      ggtitle("Histograms of estimated hrs") + 
      xlab("Hazard ratio estimates")
  }
  
  pl.obj    
}


plot.band<-function(x,mean.value,lower,upper,show.axes=F,band=TRUE,ltype="l",lty=1,xlabel=NULL,ylabel=NULL,color="grey",ylim=c(min(lower,na.rm=TRUE),max(upper,na.rm=TRUE))){
  plot(x[order(x)],mean.value[order(x)],type="n",axes=show.axes,xlab=xlabel,lty=lty,
       ylab=ylabel,ylim=ylim)
  if(band){
    polygon(c(x[order(x)],rev(x[order(x)])),c(lower[order(x)],rev(upper[order(x)])),col=color,border=FALSE)
    lines(x[order(x)],mean.value[order(x)],lty=lty,lwd=2.5,type=ltype)
  }
}


# function to prettify output
# stolen from Björn Bornkamp and Kaspar Rufibach
#https://oncoestimand.github.io/princ_strat_drug_dev/princ_strat_example.html
# Modified to work with beamer 

prettyCox <- function(mod, dig.coef = 3, dig.p = 1){
  tab1 <- data.frame(cbind(summary(mod)$conf.int), summary(mod)$coefficients)
  tab1 <- data.frame(cbind(apply(tab1[, c(5, 1)], 1:2, disp, 
                                 dig.coef), displayCI(as.matrix(tab1[, c(3, 4)]), digit = dig.coef), 
                           formatPval2(tab1[, 9], dig.p)))
  tab1 <- data.frame(cbind(rownames(summary(mod)$conf.int), tab1))
  colnames(tab1) <- c("Variable", "Coefficient", "HR", "95\\% CI", "$p$-value")
  kable(tab1, row.names = FALSE,escape=FALSE)
}
# Functions from reporttools

formatPval2<-function (pv, digits = max(1, getOption("digits") - 2), 
                       eps = 1e-04, na.form = "NA", scientific = FALSE, includeEquality = FALSE) 
{
  if ((has.na <- any(ina <- is.na(pv)))) {
    pv <- pv[!ina]
  }
  r <- character(length(is0 <- pv < eps))
  if (any(!is0)) {
    rr <- pv <- pv[!is0]
    expo <- floor(log10(ifelse(pv > 0, pv, 1e-50)))
    fixp <- expo >= -3 | (expo == -4 & digits > 1)
    if (any(fixp)) {
      rr[fixp] <- format(pv[fixp], digits = digits, scientific = scientific)
      rr[fixp] <- disp(pv[fixp], 2, 2)
    }
    if (any(!fixp)) {
      rr[!fixp] <- format(pv[!fixp], digits = digits, scientific = scientific)
      rr[!fixp] <- disp(pv[!fixp], 2, 2)
    }
    r[!is0] <- rr
  }
  if (any(is0)) {
    digits <- max(1, digits - 2)
    if (any(!is0)) {
      nc <- max(nchar(rr, type = "w"))
      if (digits > 1 && digits + 6 > nc) {
        digits <- max(1, nc - 7)
      }
    }
    r[is0] <- format(eps, digits = digits, scientific = scientific)
  }
  frontEqual <- if (includeEquality) 
    "= "
  else ""
  r <- paste(ifelse(is0, "\\ $<$", frontEqual), r, sep = "")
  if (has.na) {
    rok <- r
    r <- character(length(ina))
    r[!ina] <- rok
    r[ina] <- na.form
  }
  return(r)
}

disp<-function (n, d1 = 2, d2 = 1) 
{
  n <- as.numeric(n)
  ind.na <- is.na(n) == FALSE
  ind <- (abs(n) >= 10^-d1)
  n[ind.na & ind] <- format(round(as.numeric(n[ind.na & ind]), 
                                  d1), nsmall = d1)
  tmp <- n[(ind.na & ind) == FALSE]
  for (i in 1:length(tmp)) {
    tmp[i] <- format(as.numeric(tmp[i]), digits = d2, scientific = FALSE)
  }
  n[(ind.na & ind) == FALSE] <- tmp
  return(n)
}

displayCI<-function (ci, digit = 2, unit = "", text = "none") 
{
  ci <- format(round(ci, digit), nsmall = digit)
  d <- 1
  if (is.matrix(ci) == TRUE) {
    d <- nrow(ci)
  }
  else {
    ci <- matrix(ci, ncol = 2)
  }
  ci <- sub(" ", "", ci)
  disp.ci <- rep(NA, d)
  for (i in 1:d) {
    if (identical(text, "none")) {
      string <- paste("[", ci[i, 1], unit, ", ", 
                      ci[i, 2], unit, "]", sep = "")
    }
    if (identical(text, "german")) {
      string <- paste("von ", ci[i, 1], unit, " bis ", 
                      ci[i, 2], unit, sep = "")
    }
    if (identical(text, "english")) {
      string <- paste("from ", ci[i, 1], unit, " to ", 
                      ci[i, 2], unit, sep = "")
    }
    disp.ci[i] <- string
  }
  return(disp.ci)
}




tabcoxph2<-function(fit,
                    columns = c("beta.se", "hr.ci", "p"),
                    var.labels = NULL,
                    factor.compression = 1,
                    sep.char = ", ",
                    decimals = 2,
                    formatp.list = NULL) {
  
  # Error checking
  if (! "coxph" %in% class(fit)) {
    stop("The input 'fit' must be a fitted 'coxph'.")
  }
  if (! is.null(columns) &&
      ! all(columns %in% c("beta", "se", "betaci", "beta.se", "beta.ci", "or",
                           "hr", "hrci", "hr.ci", "z", "p"))) {
    stop("Each element of 'columns' must be one of the following: 'beta', 'se', 'betaci', 'beta.se', 'beta.ci', 'hr', 'hrci', 'hr.ci', 'z', 'p'.")
  }
  if (! factor.compression %in% 1: 5) {
    stop("The input 'factor.compression' must be set to 1, 2, 3, 4, or 5.")
  }
  if (! is.character(sep.char)) {
    stop("The input 'sep.char' must be a character string.")
  }
  if (! (is.numeric(decimals) && decimals >= 0 &&
         decimals == as.integer(decimals))) {
    stop("The input 'decimals' must be a non-negative integer.")
  }
  if (! is.null(formatp.list) &&
      ! (is.list(formatp.list) && all(names(formatp.list) %in%
                                      names(as.list(args(formatp)))))) {
    stop("The input 'format.p' must be a named list of arguments to pass to 'formatp'.")
  }
  
  # Extract info from fit
  invisible(capture.output(summary.fit <- summary(fit)))
  coefmat <- summary.fit$coefficients
  rownames.coefmat <- rownames(coefmat)
  betas <- coefmat[, "coef"]
  hrs <- coefmat[, "exp(coef)"]
  ses <- coefmat[, "se(coef)"]
  zs <- coefmat[, "z"]
  ps <- coefmat[, "Pr(>|z|)"]
  confint.fit <- confint(fit)
  lower <- confint.fit[, 1]
  upper <- confint.fit[, 2]
  
  # Convert decimals to variable for sprintf
  spf <- paste("%0.", decimals, "f", sep = "")
  
  # Initialize table
  df <- data.frame(Variable = rownames.coefmat, stringsAsFactors = FALSE)
  
  # Loop through and add columns requested
  for (column in columns) {
    
    if (column == "beta") {
      
      df$`Beta` <- sprintf(spf, betas)
      
    } else if (column == "se") {
      
      df$`SE` <- sprintf(spf, ses)
      
    } else if (column == "betaci") {
      
      df$`95% CI` <- paste("(", sprintf(spf, lower), sep.char,
                           sprintf(spf, upper), ")", sep = "")
      
    } else if (column == "beta.se") {
      
      df$`Beta (SE)` <- paste(sprintf(spf, betas), " (",
                              sprintf(spf, ses), ")", sep = "")
      
    } else if (column == "beta.ci") {
      
      df$`Beta (95% CI)` <- paste(sprintf(spf, betas), " (",
                                  sprintf(spf, lower), sep.char,
                                  sprintf(spf, upper), ")", sep = "")
      
    }  else if (column == "hr") {
      
      df$`HR` <- sprintf(spf, exp(betas))
      
    } else if (column == "hrci") {
      
      df$`95% CI` <- paste("(", sprintf(spf, exp(lower)), sep.char,
                           sprintf(spf, exp(upper)), ")", sep = "")
      
    } else if (column == "hr.ci") {
      
      df$`HR (95% CI)` <- paste(sprintf(spf, exp(betas)), " (",
                                sprintf(spf, exp(lower)), sep.char,
                                sprintf(spf, exp(upper)), ")", sep = "")
      
    } else if (column == "z") {
      
      df$`z` <- sprintf(spf, zs)
      
    } else if (column == "p") {
      
      df$`P` <- do.call(formatp, c(list(p = ps), formatp.list))
      
    }
    
  }
  
  # Clean up factor variables
  spaces <- ""
  xlevels <- fit$xlevels
  if (length(xlevels) > 0) {
    for (ii in 1: length(xlevels)) {
      varname.ii <- names(xlevels)[ii]
      levels.ii <- xlevels[[ii]]
      locs <- which(df$Variable %in% paste(varname.ii, levels.ii, sep = ""))
      if (factor.compression == 1) {
        
        # Rows are Variable, Level 1 (ref), Level 2, ...
        df$Variable[locs] <- gsub(pattern = varname.ii, replacement = spaces,
                                  x = df$Variable[locs], fixed = TRUE)
        newrows <- matrix("", nrow = 2, ncol = ncol(df), dimnames = list(NULL, names(df)))
        newrows[2, ] <- ""
        newrows[1, 1] <- ifelse(varname.ii %in% names(var.labels), var.labels[[varname.ii]], varname.ii)
        newrows[2, 1] <- paste(spaces, paste(levels.ii[1], " (ref)", sep = ""), sep = "")
        df <- rbind(df[setdiff(1: locs[1], locs[1]), ], newrows, df[locs[1]: nrow(df), ])
        
      } else if (factor.compression == 2) {
        
        # Rows are Variable (ref = Level 1), Level 2, ...
        df$Variable[locs] <- gsub(pattern = varname.ii, replacement = spaces, x = df$Variable[locs])
        newrow <- matrix("", nrow = 1, ncol = ncol(df), dimnames = list(NULL, names(df)))
        newrow[1, 1] <- paste(
          ifelse(varname.ii %in% names(var.labels), var.labels[[varname.ii]], varname.ii),
          " (ref = ", levels.ii[1], ")", sep = ""
        )
        df <- rbind(df[setdiff(1: locs[1], locs[1]), ], newrow, df[locs[1]: nrow(df), ])
        
      } else if (factor.compression == 3) {
        
        # Rows are Level 1 (ref), Level 2, ...
        df$Variable[locs] <- gsub(pattern = varname.ii, replacement = "", x = df$Variable[locs])
        newrow <- matrix("\textendash;", nrow = 1, ncol = ncol(df), dimnames = list(NULL, names(df)))
        newrow[1, 1] <- paste(levels.ii[1], "", sep = "")
        df <- rbind(df[setdiff(1: locs[1], locs[1]), ], newrow, df[locs[1]: nrow(df), ])
        
      } else if (factor.compression == 4) {
        
        # Rows are Level 2 (ref = Level 1), ...
        df$Variable[locs] <- paste(
          gsub(pattern = varname.ii, replacement = "", x = df$Variable[locs]),
          " (ref = ", levels.ii[1], ")", sep = ""
        )
        
      } else if (factor.compression == 5) {
        
        # Rows are Level 2, ...
        df$Variable[locs] <- gsub(pattern = varname.ii, replacement = "", x = df$Variable[locs])
        
      }
    }
  }
  
  # Clean up interaction terms
  interactions <- grep(":", attr(fit$terms, "term.labels"), value = TRUE)
  for (interaction.ii in interactions) {
    components <- unlist(strsplit(interaction.ii, ":"))
    locs <- intersect(grep(components[1], df$Variable),
                      grep(components[2], df$Variable))
    if (length(locs) == 1) {
      components <- c(
        ifelse(components[1] %in% names(var.labels), var.labels[[components[1]]], components[1]),
        ifelse(components[2] %in% names(var.labels), var.labels[[components[2]]], components[2])
      )
      df$Variable[locs] <- paste(components, collapse = " by ")
    } else {
      labs <- df$Variable[locs]
      labs <- gsub(components[1], "", labs)
      labs <- gsub(components[2], "", labs)
      labs <- gsub(":", ", ", labs)
      df$Variable[locs] <- paste(spaces, labs, sep = "")
      newrow <- matrix("", nrow = 1, ncol = ncol(df), dimnames = list(NULL, names(df)))
      components <- c(
        ifelse(components[1] %in% names(var.labels), var.labels[[components[1]]], components[1]),
        ifelse(components[2] %in% names(var.labels), var.labels[[components[2]]], components[2])
      )
      newrow[1, 1] <- paste(components, collapse = " by ")
      df <- rbind(df[setdiff(1: locs[1], locs[1]), ], newrow, df[locs[1]: nrow(df), ])
    }
  }
  
  # Clean up polynomial terms
  polynomials <- grep("poly(", attr(fit$terms, "term.labels"), fixed = TRUE, value = TRUE)
  for (polynomial.ii in polynomials) {
    split.ii <- unlist(strsplit(polynomial.ii, split = ", "))
    varname.ii <- substring(split.ii[1], first = 6)
    poly.order <- as.numeric(split.ii[2])
    locs <- grep(polynomial.ii, df$Variable, fixed = TRUE)
    varname.ii <- ifelse(varname.ii %in% names(var.labels), var.labels[[varname.ii]], varname.ii)
    if (poly.order == 1) {
      df$Variable[locs] <- varname.ii
    } else if (poly.order == 2) {
      df$Variable[locs] <- c(varname.ii, paste(varname.ii, "squared"))
    } else if (poly.order == 3) {
      df$Variable[locs] <- c(varname.ii, paste(varname.ii, c("squared", "cubed")))
    } else {
      df$Variable[locs] <- c(varname.ii, paste(varname.ii, 2: poly.order, sep = "^"))
    }
  }
  
  # Add user-specified labels for numeric variables
  if (! is.null(var.labels)) {
    dataClasses <- attr(fit$terms, "dataClasses")
    numerics <- names(dataClasses[dataClasses == "numeric"])
    if (length(numerics) > 0) {
      for (varname.ii in numerics) {
        loc <- which(df$Variable == varname.ii)
        if (length(loc) == 1) {
          df$Variable[loc] <- ifelse(varname.ii %in% names(var.labels),
                                     var.labels[[varname.ii]], varname.ii)
        }
      }
    }
  }
  
  # Remove row names and return table
  rownames(df) <- NULL
  return(df)
  # return(df %>% kable(escape = FALSE) %>% kable_styling(full_width = FALSE))
  
}
