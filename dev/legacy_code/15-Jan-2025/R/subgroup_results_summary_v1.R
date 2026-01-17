
hrCI_format<-function(hrest,locs=c(1,2,3),digs=2){
  # locs is the location of est, lower, upper
  hrstat<-round(hrest,digs)
  a<-paste0(hrstat[1]," (")
  a<-paste0(a,hrstat[2])
  a<-paste0(a,",")
  a<-paste0(a,hrstat[3])
  a<-paste0(a,")")
  return(c(a))
}


SG_table<-function(df,fs_bc=NULL,SG_flag,outcome.name,event.name,treat.name,strata.name=NULL,sg1_name=NULL,sg0_name=NULL,draws=0, hr_0a=NA, hr_1a=NA,details=FALSE,
                   return_medians=TRUE,est.scale="hr"){
  
if(est.scale=="1/hr"){
df$treat2 <- 1-df[,treat.name]
treat.name <- c("treat2")
}  
  
N <- nrow(df)
  
if(!is.null(fs_bc)){  
  
  resH <- fs_bc$H_estimates[,c("H2","H2_lower","H2_upper")]
  # Describing CIs
  Hstat<-c(unlist(resH[1,]))[c(1,2,3)]
  Hstat<-round(Hstat,2)
  a<-paste0(Hstat[1]," (")
  a<-paste0(a,Hstat[2])
  a<-paste0(a,", ")
  a<-paste0(a,Hstat[3])
  hr_H<-paste0(a,")")
  
  # hr_1
  resH <- fs_bc$Hc_estimates[,c("H2","H2_lower","H2_upper")]
  # Describing CIs
  Hstat<-c(unlist(resH[1,]))[c(1,2,3)]
  Hstat<-round(Hstat,2)
  a<-paste0(Hstat[1]," (")
  a<-paste0(a,Hstat[2])
  a<-paste0(a,", ")
  a<-paste0(a,Hstat[3])
  hr_Hc<-paste0(a,")")
  
if(est.scale=="hr"){
hr_1a <- hr_Hc
hr_0a <- hr_H
}
  
if(est.scale=="1/hr"){
hr_1a <- hr_H
hr_0a <- hr_Hc
}
}

  if(SG_flag!="ITT"){
    df_0<-subset(df,df[,SG_flag]==0)
    df_1<-subset(df,df[,SG_flag]==1)
    
    # df_0 --> "treatment not recommended"
    Y0<-df_0[,outcome.name]
    E0<-df_0[,event.name]
    Treat0<-df_0[,treat.name]
    Strata0 <- rep("All",length(Y0))
    
    if(!is.null(strata.name)) Strata0 <- df_0[,strata.name]
    
    # Un-adjusted CIs
    hr_0<-summary(coxph(Surv(Y0,E0)~Treat0+strata(Strata0)))$conf.int[c(1,3,4)]
    hr_0<-hrCI_format(hrest=hr_0)
    
    km_0<-summary(survfit(Surv(Y0,E0)~Treat0))
    # Medians
    meds_0<-as.numeric(round(c(km_0$table[,"median"]),1))
    
    n_0<-c(length(Y0))
    pct0 <- 100*round(n_0/N,2)
    n_0 <- paste0(n_0," (")
    n_0 <- paste0(n_0,pct0)
    n_0 <- paste0(n_0,"%)")
    
    # RMST
    # return rmst0_fit
    rmst0_fit<-KM.MeanTrunc(time=Y0,delta=E0,z=Treat0,draws=draws,details=details,get.band=ifelse(draws>0,TRUE,FALSE),plotband=FALSE)
    
    if(draws>0) drmst0<-hrCI_format(hrest=c(rmst0_fit$dhat,rmst0_fit$lower.star,rmst0_fit$upper.star))
    if(draws==0) drmst0<-hrCI_format(hrest=c(rmst0_fit$dhat,rmst0_fit$lower,rmst0_fit$upper))
    
    
    # Use resampling based CIs
    
    # Output "S0" (treatment not recommended) summaries
    #res_0<-cbind(sg0_name,n_0,length(Y0[Treat0==1]),length(Y0[Treat0==0]),meds_0[2],meds_0[1],hr_0,drmst0)
    # Revise to return individual RMST's instead of medians
    
    rmst_arms <- as.numeric(round(c(rmst0_fit$m1L,rmst0_fit$m0L),1))
    
    if(!return_medians) res_0<-cbind(sg0_name,n_0,length(Y0[Treat0==1]),length(Y0[Treat0==0]),rmst_arms[1],rmst_arms[2],drmst0,hr_0,hr_0a)
    
    if(return_medians) res_0<-cbind(sg0_name,n_0,length(Y0[Treat0==1]),length(Y0[Treat0==0]),meds_0[2],meds_0[1],drmst0,hr_0,hr_0a)
      
  
    # df_1 --> "treatment recommended"
    
    Y1<-df_1[,outcome.name]
    E1<-df_1[,event.name]
    Treat1<-df_1[,treat.name]
    
    Strata1 <- rep("All",length(Y1))
    
    if(!is.null(strata.name)) Strata1 <- df_1[,strata.name]
    
    
    hr_1<-summary(coxph(Surv(Y1,E1)~Treat1+strata(Strata1)))$conf.int[c(1,3,4)]
    hr_1<-hrCI_format(hrest=hr_1)
    
    km_1<-summary(survfit(Surv(Y1,E1)~Treat1))
    # Medians
    meds_1<-as.numeric(round(c(km_1$table[,"median"]),1))
    
    n_1<-c(length(Y1))
    
    pct1 <- 100*round(n_1/N,2)
    n_1 <- paste0(n_1," (")
    n_1 <- paste0(n_1,pct1)
    n_1 <- paste0(n_1,"%)")
    
    
    rmst1_fit<-KM.MeanTrunc(time=Y1,delta=E1,z=Treat1,draws=draws,details=details,get.band=ifelse(draws>0,TRUE,FALSE),plotband=FALSE)
    
    if(draws>0) drmst1<-hrCI_format(hrest=c(rmst1_fit$dhat,rmst1_fit$lower.star,rmst1_fit$upper.star))
    if(draws==0) drmst1<-hrCI_format(hrest=c(rmst1_fit$dhat,rmst1_fit$lower,rmst1_fit$upper))
    
    # Output "S1" (not recommend) summaries
    #res_1<-cbind(sg1_name,n_1,length(Y1[Treat1==1]),length(Y1[Treat1==0]),meds_1[2],meds_1[1],hr_1a,hr_1,drmst1)
    
    rmst_arms <- as.numeric(round(c(rmst1_fit$m1L,rmst1_fit$m0L),1))
    
    if(!return_medians) res_1<-cbind(sg1_name,n_1,length(Y1[Treat1==1]),length(Y1[Treat1==0]),rmst_arms[1],rmst_arms[2],drmst1,hr_1,hr_1a)
    if(return_medians) res_1<-cbind(sg1_name,n_1,length(Y1[Treat1==1]),length(Y1[Treat1==0]),meds_1[2],meds_1[1],drmst1,hr_1,hr_1a)
    
    res_out<-rbind(res_0,res_1)
    
    colnames(res_out)<-c("Subgroup","n","n1","n0","m1","m0","RMST","HR","HR*")
    
    return(list(res_out=res_out,rmst0_fit=rmst0_fit,rmst1_fit=rmst1_fit))
  }
  
  if(SG_flag=="ITT"){
    sg_name<-c("ITT")
    Y<-df[,outcome.name]
    E<-df[,event.name]
    Treat<-df[,treat.name]
    
    Strata <- rep("All",length(Y))
    
    if(!is.null(strata.name)) Strata <- df[,strata.name]
    
    
    hr<-summary(coxph(Surv(Y,E)~Treat+strata(Strata)))$conf.int[c(1,3,4)]
    km<-summary(survfit(Surv(Y,E)~Treat))
    # Medians
    meds<-as.numeric(round(c(km$table[,"median"]),1))
    hr<-hrCI_format(hrest=hr)
    n<-c(length(Y))
    # RMST
    # return rmst_fit
    rmst_fit<-KM.MeanTrunc(time=Y,delta=E,z=Treat,draws=draws,details=details,get.band=ifelse(draws>0,TRUE,FALSE),plotband=FALSE)
    
    if(draws>0) drmst<-hrCI_format(hrest=c(rmst_fit$dhat,rmst_fit$lower.star,rmst_fit$upper.star))
    if(draws==0) drmst<-hrCI_format(hrest=c(rmst_fit$dhat,rmst_fit$lower,rmst_fit$upper))
    
    # Use resampling based CIs
    #res<-cbind(sg_name,n,length(Y[Treat==1]),length(Y[Treat==0]),meds[2],meds[1],hr,drmst)
    
    rmst_arms <- as.numeric(round(c(rmst_fit$m1L,rmst_fit$m0L),1))
    
    if(!return_medians) res <-cbind(sg_name,n,length(Y[Treat==1]),length(Y[Treat==0]),rmst_arms[1],rmst_arms[2],drmst,hr,NA)
    if(return_medians) res <-cbind(sg_name,n,length(Y[Treat==1]),length(Y[Treat==0]),meds[2],meds[1],drmst,hr,NA)
    
    
    colnames(res)<-c("Subgroup","n","n1","n0","m1","m0","RMST","HR","HR*")
    
    return(list(res_out=res,rmst_fit=rmst_fit))
  }
}



#plotband_survdiff<-function(res,xlab="Months"){
#  plot.band(x=res$dpoints,mean=res$diff.dpoints,xlabel=xlab,
#            ylabel=expression(delta(t)==hat(S)[1](t)-hat(S)[0](t)),
#            band=TRUE,lower=res$sb.lower,upper=res$sb.upper,show.axes=TRUE,ltype="l")
#  lines(res$dpoints,res$pw.upper,type="l",lty=2,lwd=0.50)
#  lines(res$dpoints,res$pw.lower,type="l",lty=2,lwd=0.50)
#  abline(h=0.0,lty=1,col="black",lwd=2)
#  rug(res$dpoints)
#}

