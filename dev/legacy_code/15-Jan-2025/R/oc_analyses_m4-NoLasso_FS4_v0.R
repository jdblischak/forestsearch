
# Default parameters
oc_analyses_m4NoLasso_FS4<-function(sim){
  
  ans.analyses<-NULL  
  
  x<-sim_aftm4_gbsg(dgm=dgm,n=N,maxFollow=maxFollow,muC.adj=muC.adj,simid=sim)

  if(n_add_noise==0){
  confounders.name <- c("z1","z2","z3","z4","z5","size","grade3")
  }
    
  if(n_add_noise==5){
  set.seed(8316951+1000*sim)
  # Add 5 noise 
  x$noise1<-rnorm(N)
  x$noise2<-rnorm(N)
  x$noise3<-rnorm(N)
  x$noise4<-rnorm(N)
  x$noise5<-rnorm(N)
  confounders.name <- c("z1","z2","z3","z4","z5","size","grade3","noise1","noise2","noise3","noise4","noise5")
  }
  
  if(n_add_noise==3){
    set.seed(8316951+1000*sim)
    # Add 3 noise 
    x$noise1<-rnorm(N)
    x$noise2<-rnorm(N)
    x$noise3<-rnorm(N)
    confounders.name <- c("z1","z2","z3","z4","z5","size","grade3","noise1","noise2","noise3")
  }
  
  
  # Initiate with # and prop H
  sizeH_true<-with(x,sum(flag.harm))
  propH_true<-with(x,mean(flag.harm))
  sizeHc_true<-with(x,sum(!flag.harm))
  propHc_true<-with(x,mean(!flag.harm))
  
  dfHpop<-data.table(sim,sizeH_true,propH_true,sizeHc_true,propHc_true)
  
  
showsim <- FALSE
# show first and last simulation
if(sim <= 1 | sim>=(Nsims-1)) showsim <- TRUE

  if(get.FS){
    use_lasso <- FALSE
    use_grf <- FALSE
    
    ansFS<-NULL
    analysis<-label.analyses[1]
   
      y.fs<-try(forestsearch(df.analysis=x, Allconfounders.name=confounders.name,
                        details=showsim,use_lasso=use_lasso, use_grf=use_grf,
                        conf_force=NULL, 
                        outcome.name=outcome.name,treat.name=treat.name,event.name=event.name,id.name=id.name,
                        n.min=nmin.fs,hr.threshold=hr.threshold,hr.consistency=hr.consistency,fs.splits=fs.splits,
                        d0.min=d.min,d1.min=d.min,pstop_futile=pstop_futile,
                        pconsistency.threshold=pconsistency.threshold, stop.threshold=stop.threshold,
                        max.minutes=max.minutes,maxk=maxk,by.risk=12,
                        plot.sg=showsim,vi.grf.min=vi.grf.min),TRUE)
    
    if(!inherits(y.fs,"try-error")){
      # sg_focus="hr" (default FS)

      ansFS<-fs.estimates.out(df=x,fs.res=y.fs$grp.consistency$out_hr,dgm=dgm,
                              cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,
                              analysis=analysis)
      ansFS<-cbind(dfHpop,ansFS)
      ans.analyses<-rbind(ans.analyses,ansFS)
    }
    
    if(inherits(y.fs,"try-error")){
      ansFS<-fs.estimates.out(df=x,fs.res=NULL,dgm=dgm,
                              cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,
                              analysis=analysis)
      ansFS<-cbind(dfHpop,ansFS)
      ans.analyses<-rbind(ans.analyses,ansFS)
      if(round(c(ansFS$hr.Hc.hat-ansFS$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (FS)")
    }
      rm("y.fs")
  }
  
  if(get.FS){
    ansFS<-NULL
    analysis<-"FSlg"
    use_lasso <- FALSE
    use_grf <- TRUE
  
    y.fs<-try(forestsearch(df.analysis=x, Allconfounders.name=confounders.name,
                           details=showsim,use_lasso=use_lasso, use_grf=use_grf,
                           conf_force=NULL, 
                           frac.tau=frac.tau_fs,
                           dmin.grf=dmin.grf_fs,
                           grf_depth=maxdepth_fs,
                           outcome.name=outcome.name,treat.name=treat.name,event.name=event.name,id.name=id.name,
                           n.min=nmin.fs,hr.threshold=hr.threshold,hr.consistency=hr.consistency,fs.splits=fs.splits,
                           d0.min=d.min,d1.min=d.min,pstop_futile=pstop_futile,
                           pconsistency.threshold=pconsistency.threshold, stop.threshold=stop.threshold,
                           max.minutes=max.minutes,maxk=maxk,by.risk=12,
                           plot.sg=showsim,vi.grf.min=vi.grf.min),TRUE)
    
    
    
   if(!inherits(y.fs,"try-error")){
      ansFS<-fs.estimates.out(df=x,fs.res=y.fs$grp.consistency$out_hr,dgm=dgm,
                              cox.formula.sim=cox.formula.sim,frac.tau_fs=frac.tau_fs,
                              cox.formula.adj.sim=cox.formula.adj.sim,
                              analysis=analysis)
      ansFS<-cbind(dfHpop,ansFS)
      ans.analyses<-rbind(ans.analyses,ansFS)
    }
    
    if(inherits(y.fs,"try-error")){
      ansFS<-fs.estimates.out(df=x,fs.res=NULL,dgm=dgm,
                              cox.formula.sim=cox.formula.sim,frac.tau_fs=frac.tau_fs,
                              cox.formula.adj.sim=cox.formula.adj.sim,
                              analysis=analysis)
      ansFS<-cbind(dfHpop,ansFS)
      ans.analyses<-rbind(ans.analyses,ansFS)
      if(round(c(ansFS$hr.Hc.hat-ansFS$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (FS)")
    }
  }
  
  # frac.tau<-1.0
  if(get.GRF){
    
    ansGRF<-NULL
    analysis<-label.analyses[2]
    
    grf.est<-try(grf.subg.harm.survival(data=x,confounders.name=confounders.name,outcome.name=outcome.name,
                                        event.name=event.name,id.name=id.name,treat.name=treat.name,n.min=n.min,
                                        dmin.grf=dmin.grf,frac.tau=1.0, maxdepth=maxdepth, details=FALSE),TRUE)
    
    if(!inherits(grf.est,"try-error")){
      ansGRF<-grf.estimates.out(df=x,grf.est=grf.est,dgm=dgm,cox.formula.sim=cox.formula.sim,
                                cox.formula.adj.sim=cox.formula.adj.sim,analysis=analysis,frac.tau=1.0)
      ansGRF<-cbind(dfHpop,ansGRF)
      ans.analyses<-rbind(ans.analyses,ansGRF)
    }
    
    if(inherits(grf.est,"try-error")){
      ansGRF<-grf.estimates.out(df=x,grf.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                                cox.formula.adj.sim=cox.formula.adj.sim,analysis=analysis,frac.tau=1.0)
      ansGRF<-cbind(dfHpop,ansGRF)
      ans.analyses<-rbind(ans.analyses,ansGRF)
      if(round(c(ansGRF$hr.Hc.hat-ansGRF$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (GRF)")
    }
    rm("grf.est")
    
  } # GRF analysis frac.tau=100%
  
  
  if(get.GRF){
    ansGRF<-NULL
    analysis<-label.analyses[7]
    
    grf.est<-try(grf.subg.harm.survival(data=x,confounders.name=confounders.name,outcome.name=outcome.name,
                                        event.name=event.name,id.name=id.name,treat.name=treat.name,n.min=n.min,
                                        dmin.grf=dmin.grf,frac.tau=frac.tau, maxdepth=maxdepth, details=FALSE),TRUE)
    
    if(!inherits(grf.est,"try-error")){
      
      ansGRF<-grf.estimates.out(df=x,grf.est=grf.est,dgm=dgm,cox.formula.sim=cox.formula.sim,
                                cox.formula.adj.sim=cox.formula.adj.sim,analysis=analysis,frac.tau=frac.tau)
      
      ansGRF<-cbind(dfHpop,ansGRF)
      
      ans.analyses<-rbind(ans.analyses,ansGRF)
    }
    
    if(inherits(grf.est,"try-error")){
      ansGRF<-grf.estimates.out(df=x,grf.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                                cox.formula.adj.sim=cox.formula.adj.sim,analysis=analysis,frac.tau=1.0)
      ansGRF<-cbind(dfHpop,ansGRF)
      ans.analyses<-rbind(ans.analyses,ansGRF)
      if(round(c(ansGRF$hr.Hc.hat-ansGRF$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (GRF)")
    }
    rm("grf.est")
  } # GRF analysis frac.tau specified
  
  
  
  if(get.VT){
    # VT
    x2<-get.FG(x,tte.name="y.sim",event.name="event.sim")
    # VT at tau1, censored version 
    tau<-24
    analysis<-label.analyses[3]
    ansVT<-NULL
    vt.estA<-try(vt.subg.harm.survival(data=x2,cf.names=confounders.name,tte.name="Y.FG",event.name="event.YFG",treat.name="treat",ntree=ntree,
                                       seed=8316951,data.pred=NULL,details=FALSE,n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau),TRUE)
    
    if(!inherits(vt.estA,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=vt.estA,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
    }
    
    if(inherits(vt.estA,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
      if(round(c(ansVT$hr.Hc.hat-ansVT$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (VT 24)")
    }
    
    
    ansVT<-NULL
    analysis<-label.analyses[4]
    # VT at tau1, NON-censored version t.sim (Note: event.YFG are 'dummy' event indicators == 1 for all)
    vt.estB<-try(vt.subg.harm.survival(data=x2,cf.names=confounders.name,tte.name="t.sim",event.name="event.YFG",treat.name="treat",ntree=ntree,
                                       seed=8316951,data.pred=NULL,details=FALSE,n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau),TRUE)
    
    if(!inherits(vt.estB,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=vt.estB,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
    }
    
    if(inherits(vt.estB,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
      if(round(c(ansVT$hr.Hc.hat-ansVT$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (VT 24 ideal)")
    }
    
    tau<-36
    ansVT<-NULL
    analysis<-label.analyses[5]
    vt.estC<-try(vt.subg.harm.survival(data=x2,cf.names=confounders.name,tte.name="Y.FG",event.name="event.YFG",treat.name="treat",ntree=ntree,
                                       seed=8316951,data.pred=NULL,details=FALSE,n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau),TRUE)
    
    if(!inherits(vt.estC,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=vt.estC,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      ans.analyses<-rbind(ans.analyses,ansVT)
    }
    
    if(inherits(vt.estC,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
      if(round(c(ansVT$hr.Hc.hat-ansVT$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (VT 36)")
    }
    
    
    ansVT<-NULL
    analysis<-label.analyses[6]
    # VT at tau1, NON-censored version t.sim (Note: event.YFG are 'dummy' event indicators == 1 for all)
    vt.estD<-try(vt.subg.harm.survival(data=x2,cf.names=confounders.name,tte.name="t.sim",event.name="event.YFG",treat.name="treat",ntree=ntree,
                                       seed=8316951,data.pred=NULL,details=FALSE,n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau),TRUE)
    if(!inherits(vt.estD,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=vt.estD,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
    }
    
    if(inherits(vt.estD,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
      if(round(c(ansVT$hr.Hc.hat-ansVT$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (VT 36 ideal)")
    }
    
  }

#t.now<-proc.time()[3]
#t.min<-round((t.now-t.start.all)/60,4)
#cat("Simulation done for sim, minutes so far=",c(sim),"\n")
return(ans.analyses)
}



tab_tests<-function(res,analyses=unique(res$analysis)){
  test<-get.summary.VTFS(result=subset(res,analysis==analyses[1]))
  names(test)<-analyses[1]
  test.all<-test
  for(aa in 2:length(analyses)){
    test<-get.summary.VTFS(result=subset(res,analysis==analyses[aa]))
    names(test)<-analyses[aa]
    test.all<-cbind(test.all,test)
  }  
  return(test.all)  
}


out.results<-function(res,dgm,output.file=NULL,t.min,details=TRUE,out_analysis="FS"){
  if(details){
    print(output.file)
    print(head(res))
    cat("Subgroup HRs: H, H^c, Causal=",c(dgm$hr.H.true,dgm$hr.Hc.true,dgm$hr.causal),"\n")
    out.analysis<-subset(res,analysis==out_analysis)
    cat("Simulations=",c(nrow(out.analysis)),"\n")
    pC<-mean(out.analysis$p.cens)
    cat("Avg censoring=",c(pC),"\n")
    cat("Min,Max,Avg tau.max=",c(min(out.analysis$taumax),max(out.analysis$taumax),mean(out.analysis$taumax)),"\n")
  }
  t.hours<-t.min/60
  
  if(!is.na(dgm$hr.H.true)){
    hrH.plim<-mean(res$hr.H.true)
    hrH.true<-dgm$hr.H.true
    # calculate at n.min=60 and Expected N (n.sg_super)
    n.sg_super<-round(mean(dgm$df.super_rand$flag.harm)*N,0)
    pAnyH.approx1 <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrH.true,pC=pC,n.sg=60,k1=log(hr.threshold),k2=log(hr.consistency))$integral
    pAnyH.approx2 <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrH.true,pC=pC,n.sg=n.sg_super,k1=log(hr.threshold),k2=log(hr.consistency))$integral
    pAnyH.approx3 <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrH.plim,pC=pC,n.sg=n.sg_super,k1=log(hr.threshold),k2=log(hr.consistency))$integral
    if(details){
      cat("P(H) approximation at causal(H), n(sg)=60, approx=",c(hrH.true, 60, pAnyH.approx1),"\n")
      cat("P(H) approximation at causal(H), Avg(n(sg)), approx=",c(hrH.true, n.sg_super, pAnyH.approx2),"\n")
      cat("P(H) approximation at plim(H), Avg(n(sg)), approx=",c(hrH.plim, n.sg_super, pAnyH.approx3),"\n")
    }
  }
  
  if(is.na(dgm$hr.H.true)){
    hrHc.plim<-mean(res$hr.Hc.true)
    hrHc.true<-dgm$hr.Hc.true
    pAnyH.approx1 <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrHc.true,pC=pC,n.sg=60,k1=log(hr.threshold),k2=log(hr.consistency))$integral
    pAnyH.approx2 <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrHc.plim,pC=pC,n.sg=60,k1=log(hr.threshold),k2=log(hr.consistency))$integral
    pAnyH.approx3<-NULL
    if(details){
      cat("P(H) approximation at causal(Hrc), n=60, approx=",c(hrHc.true,pAnyH.approx1),"\n")
      cat("P(H) approximation at plim(Hrc), n=60, approx=",c(hrHc.plim,pAnyH.approx2),"\n")
    }
  }
  
  
  # Output file
  oc.tab<-tab_tests(res=res)
  
  if(details){
    cat("Minutes,hours",c(t.min,t.hours),"\n")
    print(oc.tab)
  }
  
  if(!is.null(output.file)) save(pAnyH.approx1,pAnyH.approx2,pAnyH.approx3,t.hours,dgm,res,oc.tab,file=output.file)
  
  return(list(pAnyH.approx1=pAnyH.approx1,pAnyH.approx2=pAnyH.approx2,pAnyH.approx3=pAnyH.approx3,  
              t.hours,oc.tab=oc.tab))
}

get.summary.VTFS<-function(result,analysis="VT",dig1=2,dig2=3){
  # summarize finding H
  # remove "found.al3" from table
  dfH1a<-result[,c("any.H","ppv","npv","sensitivity","specificity")]  
  temp1a<-as.numeric(apply(dfH1a,2,mean,na.rm=TRUE))
  temp1a<-round(temp1a,digits=dig1)
  
  # Where hat(H) exists --> hat(Hc) exists
  resH<-subset(result,any.H==1)
  # resH includes all sequences where a subgroup is found
  
  #dfH1b<-as.matrix(resH[,c("npv")])  
  #temp1b<-as.numeric(apply(dfH1b,2,mean,na.rm=TRUE))
  #temp1b<-round(temp1b,digits=dig1)
  # removing ppv(H|H found)
  temp1<-c(temp1a)
  # Size of H
  avg.H<-mean(resH$size.H,na.rm=TRUE)
  min.H<-min(resH$size.H,na.rm=TRUE)
  max.H<-max(resH$size.H,na.rm=TRUE)
  # Hc
  avg.Hc<-mean(result$size.Hc,na.rm=TRUE)
  min.Hc<-min(result$size.Hc,na.rm=TRUE)
  max.Hc<-max(result$size.Hc,na.rm=TRUE)
  
  temp2<-round(c(avg.H,min.H,max.H,avg.Hc,min.Hc,max.Hc),digits=0)
  
  #dfH2<-result[,c("hr.H.true","hr.H.hat")]
  dfH2<-resH[,c("hr.H.true","hr.H.hat")]
  #dfH2<-subset(dfH2,any.H==1)
  temp3<-as.numeric(apply(dfH2,2,mean,na.rm=TRUE))
  temp3<-round(temp3,digits=dig2)
  
  dfH3<-resH[,c("hr.Hc.true","hr.Hc.hat")]  
  temp4<-as.numeric(apply(dfH3,2,mean,na.rm=TRUE))
  temp4<-round(temp4,digits=dig2)
  
  # Include hr.H.true and hr.Hc.true for all samples
  dfH4<-result[,c("hr.H.true","hr.Hc.true","hr.itt","hr.adj.itt")]  
  temp5<-as.numeric(apply(dfH4,2,mean,na.rm=TRUE))
  temp5<-round(temp5,digits=dig2)
  
  res.out<-c(temp1,temp2,temp3,temp4,temp5)
  # remove "found.al3"
  res.out2<-data.frame(res.out)
  
  # Rename ppv=ppH, npv=ppHc, sens=sensH, spec=sensHc
  # Swap defition of ppv and sens
  # Switch ppH <-> sensH
  # Remove hathH_ppHc
  row.names(res.out2)<-c("any.H","sensH","sensHc","ppH","ppHc",
                         "Avg(#H)","minH","maxH","Avg(#Hc)","minHc","maxHc",
                         "hat(H*)","hat(hat[H])","hat(Hc*)","hat(hat[Hc])","hat(H*)all","hat(Hc*)all","hat(ITT)all","hat(ITTadj)all")
  return(res.out2)
}



get_tabsim<-function(missC,get_stats=c("any.H","sensH","sensHc","ppH","ppHc","Avg(#H)","minH","maxH","Avg(#Hc)","minHc","maxHc"), pA,est_names,stat_names,mod.harm,Nsims){
  #pA1<-paste(pA,"/%")
  #eval(parse(text=pA))
  missC<-missC[c(get_stats),]
  missC_alt<-missC
  colnames(missC_alt)<-est_names
  rownames(missC_alt)<-stat_names
  missC<-format(missC_alt,digits=3,drop0trailing = TRUE)
  # For null model, set sensitivity to missing
  if(mod.harm=="null") missC[2,]<-NA
  options(knitr.kable.NA = '.',format="latex")
  tabsim_missC<-kbl(missC,longtable=FALSE,align='c',format="latex", digits=3, booktabs=TRUE,escape=F,
                    caption="Average classification rates: ${avg}\\vert \\hat{H} \\vert$, ${min}\\vert \\hat{H} \\vert$, and ${max}\\vert \\hat{H} \\vert$, denote the average, 
minimum, and maximum of the number of subjects in the estimated subgroup $\\hat{H}$ (analogously for $\\hat{H}^{c}$).  
Note that under the null $sens(\\hat{H})$ is undefined and $ppv(\\hat{H})=0$.") %>%
    kable_styling("striped", full_width = F, latex_options="hold_position", font_size=9) %>%
    group_rows("Finding H", 1, 5) %>%
    group_rows("Size of H and H-complement", 6, 11) %>%
    footnote(general=paste('Probability approximation=', pA, '.'), footnote_as_chunk=T) %>%
    footnote(general=paste('Number of simulations=', Nsims, '.'), footnote_as_chunk=T) 
  return(tabsim_missC)
}



hr_threshold_FS<-function(n.sg,pC=0.3,theta,hr.threshold,details=FALSE,alpha=0.025){
  d.sg<-n.sg*(1-pC)
  result<-1-pnorm(sqrt(d.sg/4)*(log(hr.threshold)-log(theta)-(2/sqrt(d.sg))*qnorm(1-alpha)))
  
  if(details) cat("Approximation=",c(result),"\n")
  return(result)
}



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




