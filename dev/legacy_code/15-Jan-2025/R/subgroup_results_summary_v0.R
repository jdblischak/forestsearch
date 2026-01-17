
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


n_pcnt <- function(x,denom){
  m <-c(length(x))
  pct <- 100*round(m/denom,2)
  m <- paste0(m," (")
  m <- paste0(m,pct)
  m <- paste0(m,"%)")
  return(m)
}



SG_tab_estimates<-function(df,SG_flag,outcome.name,event.name,treat.name,strata.name=NULL,
                           sg1_name=NULL,sg0_name=NULL,draws=0, 
                           details=FALSE,return_medians=TRUE,est.scale="hr"){
  
  if(est.scale=="1/hr"){
    df$treat2 <- 1-df[,treat.name]
    treat.name <- c("treat2")
  }  
  
  N <- nrow(df)
  
  if(SG_flag!="ITT"){
    
    if(est.scale=="1/hr"){
      df_0<-subset(df,df[,SG_flag]==1)
      df_1<-subset(df,df[,SG_flag]==0)
    }  
    
    if(est.scale=="hr"){
      df_0<-subset(df,df[,SG_flag]==0)
      df_1<-subset(df,df[,SG_flag]==1)
    }  
    
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
    
    # n_0 is number in "not recommended"
    
    n_0 <- n_pcnt(x=Y0,denom=N)
    
    # RMST
    # return rmst0_fit
    rmst0_fit<-KM.MeanTrunc(time=Y0,delta=E0,z=Treat0,draws=draws,details=details,get.band=ifelse(draws>0,TRUE,FALSE),plotband=FALSE)
    
    if(draws>0) drmst0<-hrCI_format(hrest=c(rmst0_fit$dhat,rmst0_fit$lower.star,rmst0_fit$upper.star))
    if(draws==0) drmst0<-hrCI_format(hrest=c(rmst0_fit$dhat,rmst0_fit$lower,rmst0_fit$upper))
    
    n_0_treat <- n_pcnt(x=Y0[Treat0==1],denom=length(Y0))
    
    # Events
    
    d_0 <- n_pcnt(x=Y0[E0==1],denom=length(Y0))
    
    # Use resampling based CIs
    
    # Output "S0" (treatment not recommended) summaries
    #res_0<-cbind(sg0_name,n_0,length(Y0[Treat0==1]),length(Y0[Treat0==0]),meds_0[2],meds_0[1],hr_0,drmst0)
    # Revise to return individual RMST's instead of medians
    
    rmst_arms <- as.numeric(round(c(rmst0_fit$m1L,rmst0_fit$m0L),1))
    
    if(!return_medians) res_0<-cbind(sg0_name,n_0,n_0_treat,d_0,rmst_arms[1],rmst_arms[2],drmst0,hr_0)
    
    if(return_medians) res_0<-cbind(sg0_name,n_0,n_0_treat,d_0,meds_0[2],meds_0[1],drmst0,hr_0)
    
    
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
    
    
    n_1 <- n_pcnt(x=Y1,denom=N)
    
    rmst1_fit<-KM.MeanTrunc(time=Y1,delta=E1,z=Treat1,draws=draws,details=details,get.band=ifelse(draws>0,TRUE,FALSE),plotband=FALSE)
    
    if(draws>0) drmst1<-hrCI_format(hrest=c(rmst1_fit$dhat,rmst1_fit$lower.star,rmst1_fit$upper.star))
    if(draws==0) drmst1<-hrCI_format(hrest=c(rmst1_fit$dhat,rmst1_fit$lower,rmst1_fit$upper))
    
    n_1_treat <- n_pcnt(x=Y1[Treat1==1],denom=length(Y1))
    
    # Events
    
    d_1 <- n_pcnt(x=Y1[E1==1],denom=length(Y1))
    
    # Output "S1" (not recommend) summaries
    #res_1<-cbind(sg1_name,n_1,length(Y1[Treat1==1]),length(Y1[Treat1==0]),meds_1[2],meds_1[1],hr_1a,hr_1,drmst1)
    
    rmst_arms <- as.numeric(round(c(rmst1_fit$m1L,rmst1_fit$m0L),1))
    
    if(!return_medians)  res_1<-cbind(sg1_name,n_1,n_1_treat,d_1,rmst_arms[1],rmst_arms[2],drmst1,hr_1)
    if(return_medians) res_1<-cbind(sg1_name,n_1,n_1_treat,d_1,meds_1[2],meds_1[1],drmst1,hr_1)
    
    res <-rbind(res_0,res_1)
    
    colnames(res)<-c("Subgroup","n","n1","events","m1","m0","RMST","HR (95% CI)")
    
    return(res)
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
    
    n_itt <- n_pcnt(x=Y,denom=n)
    
    # RMST
    # return rmst_fit
    rmst_fit<-KM.MeanTrunc(time=Y,delta=E,z=Treat,draws=draws,details=details,get.band=ifelse(draws>0,TRUE,FALSE),plotband=FALSE)
    
    if(draws>0) drmst<-hrCI_format(hrest=c(rmst_fit$dhat,rmst_fit$lower.star,rmst_fit$upper.star))
    if(draws==0) drmst<-hrCI_format(hrest=c(rmst_fit$dhat,rmst_fit$lower,rmst_fit$upper))
    
    # Use resampling based CIs
    #res<-cbind(sg_name,n,length(Y[Treat==1]),length(Y[Treat==0]),meds[2],meds[1],hr,drmst)
    
    rmst_arms <- as.numeric(round(c(rmst_fit$m1L,rmst_fit$m0L),1))
    
    n_treat <- n_pcnt(x=Y[Treat==1],denom=n)
    d <- n_pcnt(x=Y[E==1],denom=n)
    
    if(!return_medians) res <-cbind(sg_name,n_itt,n_treat,d,rmst_arms[1],rmst_arms[2],drmst,hr)
    if(return_medians) res <-cbind(sg_name,n_itt,n_treat,d,meds[2],meds[1],drmst,hr)
    
    
    colnames(res)<-c("Subgroup","n","n1","events","m1","m0","RMST","HR (95% CI)")
    
    return(res)
  }
}


sg_tables <- function(fs, which_df="est", est_caption="Training data estimates"){
  if(which_df=="est") df <- fs$df.est
  if(which_df=="testing") df <- fs$df.test
  # ITT estimates  
  temp <- SG_tab_estimates(df=df, est.scale=fs$est.scale, SG_flag="ITT",draws=0, details=FALSE,
                           outcome.name=fs$outcome.name,event.name=fs$event.name,treat.name=fs$treat.name)
  aa <- as.data.frame(temp)
  
  if(fs$est.scale=="hr"){
    temp <- SG_tab_estimates(df=df,SG_flag="treat.recommend",draws=0, details=FALSE,
                             outcome.name=fs$outcome.name,event.name=fs$event.name,treat.name=fs$treat.name,
                             est.scale="hr",sg0_name="Questionable",sg1_name="Recommend")
    bb <- as.data.frame(temp)
  }
  
  if(fs$est.scale=="1/hr"){
    #df[,"treat.recommend2"] <- 1 -df[,"treat.recommend"]
    # treat.recommend2: 1 -- > recommend
    temp <- SG_tab_estimates(df=df,SG_flag="treat.recommend",draws=0, details=FALSE,
                             outcome.name=fs$outcome.name,event.name=fs$event.name,treat.name=fs$treat.name,
                             est.scale="1/hr",sg0_name="Questionable",sg1_name="Recommend")
    bb <- as.data.frame(temp)
  }
  
  tab_estimates <- rbind(aa,bb)
  
  tab_estimates <- gt(tab_estimates, caption=est_caption)
  
  if(fs$sg_focus=="hr") sg10  <- as.data.frame(fs$grp.consistency$out_hr$result)
  if(fs$sg_focus %in% c("minSG","hrMinSG")) sg10  <- as.data.frame(fs$grp.consistency$out_minSG$result)
  if(fs$sg_focus %in% c("maxSG","hrMaxSG")) sg10  <- as.data.frame(fs$grp.consistency$out_maxSG$result)
  
  if(fs$maxk==1){
    # Subgroup identification ="M.1"
    # Consistency rate = "Pcons"
    # hr = "hr"
    # sample size = "N"
    # Total events = "E"
    sg10 <- sg10[,c("M.1","N","E","hr","Pcons")]
    
    if(fs$est.scale=="1/hr"){
      sg10$hr <- c(1/sg10$hr)  
    }
    
    sg10_out <- sg10 |> gt() |>
      fmt_number(columns=c(4,5),decimals=6) |>
      fmt_number(columns=c(2,3),decimals=0) |>
      tab_header(title="Subgroups formed by single-factors",
                 subtitle="maxk=1")
  }
  if(fs$maxk==2){
    sg10 <- sg10[,c("M.1","M.2","N","E","hr","Pcons")]
    
    if(fs$est.scale=="1/hr"){
      sg10$hr <- c(1/sg10$hr)  
    }
    
    sg10_out <- sg10 |> gt() |>
      fmt_number(columns=c(5,6),decimals=6) |>
      fmt_number(columns=c(2,3),decimals=0) |>
      tab_header(title="Subgroups formed by two-factors",
                 subtitle="maxk=2")
  }
  
  
  return(list(tab_estimates=tab_estimates,sg10_out=sg10_out))
}



# Plotting ITT and identified subgroups

plot_found.subgroups <- function(fs, df,showITT=TRUE, byrisk=6, xlab="Months", ylab="Overall Survival", 
                                 arms=c("Experimental","Control"),title_itt=c("AP ITT"),title_benefit=c("'Consistent Benefit'")){
  
  # Note: add check for names in dataset df  
  sg.harm <- fs$sg.harm
  outcome.name <- fs$outcome.name
  treat.name <- fs$treat.name
  event.name <- fs$event.name
  est.scale <- fs$est.scale  
  
  if(est.scale=="1/hr"){
    df$treat2 <- 1-df[,treat.name]
    treat.name <- c("treat2")
  }  
  
  id_harm <- paste(sg.harm,collapse=" & ")
  
  if(showITT){
    layout.matrix <- matrix(c(1,1,2,3), 2,2, byrow=TRUE)
    layout(mat=layout.matrix, heights=c(2.25,2), widths=c(3,3))
    
    kmH.fit<-KM.plot.2sample.weighted(df=df, 
                                      tte.name=outcome.name, event.name=event.name, treat.name=treat.name,
                                      risk.set=TRUE, by.risk=byrisk,
                                      Xlab=xlab,Ylab=ylab,details=FALSE,
                                      put.legend.lr="left", put.legend.arms="top",
                                      xmed.offset=0,
                                      arms=arms)
    title(main=title_itt)
    
    # "Identified (H) subgroup" flagged via treat.recommend==0 
    if(est.scale=="hr"){
      kmH.fit<-KM.plot.2sample.weighted(df=subset(df,treat.recommend==0), 
                                        tte.name=outcome.name, event.name=event.name, treat.name=treat.name,
                                        risk.set=TRUE, by.risk=byrisk,
                                        Xlab=xlab,Ylab=ylab,details=FALSE,
                                        show.logrank=FALSE, 
                                        xmed.offset=3,
                                        put.legend.arms="top", 
                                        arms=arms)
      title(main=c(id_harm))
      
      # Strongly (consistently) positive treat.recommend==1
      kmH.fit<-KM.plot.2sample.weighted(df=subset(df,treat.recommend==1), 
                                        tte.name=outcome.name, event.name=event.name, treat.name=treat.name,
                                        risk.set=TRUE, by.risk=byrisk,
                                        Xlab=xlab,Ylab=ylab,details=FALSE, 
                                        show.logrank=FALSE, 
                                        xmed.offset=3,
                                        put.legend.arms="top", 
                                        arms=arms)
      title(main=title_benefit)
    }
    
    if(est.scale=="1/hr"){
      kmH.fit<-KM.plot.2sample.weighted(df=subset(df,treat.recommend==1), 
                                        tte.name=outcome.name, event.name=event.name, treat.name=treat.name,
                                        risk.set=TRUE, by.risk=byrisk,
                                        Xlab=xlab,Ylab=ylab,details=FALSE,
                                        show.logrank=FALSE, 
                                        xmed.offset=3,
                                        put.legend.arms="top", 
                                        arms=arms)
      title(main=c("Questionable"))
      
      # Strongly (consistently) positive treat.recommend==1
      kmH.fit<-KM.plot.2sample.weighted(df=subset(df,treat.recommend==0), 
                                        tte.name=outcome.name, event.name=event.name, treat.name=treat.name,
                                        risk.set=TRUE, by.risk=byrisk,
                                        Xlab=xlab,Ylab=ylab,details=FALSE, 
                                        show.logrank=FALSE, 
                                        xmed.offset=3,
                                        put.legend.arms="top", 
                                        arms=arms)
      title(main=title_benefit,sub=c(id_harm))
    }
  }
  
  if(!showITT){
    par(mfrow=c(1,2))
    if(est.scale=="hr"){  
      # "Identified (H) subgroup" flagged via treat.recommend==0 
      kmH.fit<-KM.plot.2sample.weighted(df=subset(df,treat.recommend==0), 
                                        tte.name=outcome.name, event.name=event.name, treat.name=treat.name,
                                        risk.set=TRUE, by.risk=byrisk,
                                        Xlab=xlab,Ylab=ylab,details=FALSE,
                                        show.logrank=FALSE, 
                                        xmed.offset=3,
                                        put.legend.arms="top", 
                                        arms=arms)
      title(main=c(id_harm))
      # Strongly (consistently) positive treat.recommend==1
      kmH.fit<-KM.plot.2sample.weighted(df=subset(df,treat.recommend==1), 
                                        tte.name=outcome.name, event.name=event.name, treat.name=treat.name,
                                        risk.set=TRUE, by.risk=byrisk,
                                        Xlab=xlab,Ylab=ylab,details=FALSE, 
                                        show.logrank=FALSE, 
                                        xmed.offset=3,
                                        put.legend.arms="top", 
                                        arms=arms)
      title(main=title_benefit)
    }
    if(est.scale=="1/hr"){  
      # "Identified (H) subgroup" flagged via treat.recommend==0 
      kmH.fit<-KM.plot.2sample.weighted(df=subset(df,treat.recommend==1), 
                                        tte.name=outcome.name, event.name=event.name, treat.name=treat.name,
                                        risk.set=TRUE, by.risk=byrisk,
                                        Xlab=xlab,Ylab=ylab,details=FALSE,
                                        show.logrank=FALSE, 
                                        xmed.offset=3,
                                        put.legend.arms="top", 
                                        arms=arms)
      title(main=c("Questionable"))
      # Strongly (consistently) positive treat.recommend==1
      kmH.fit<-KM.plot.2sample.weighted(df=subset(df,treat.recommend==0), 
                                        tte.name=outcome.name, event.name=event.name, treat.name=treat.name,
                                        risk.set=TRUE, by.risk=byrisk,
                                        Xlab=xlab,Ylab=ylab,details=FALSE, 
                                        show.logrank=FALSE, 
                                        xmed.offset=3,
                                        put.legend.arms="top", 
                                        arms=arms)
      title(main=title_benefit,sub=c(id_harm))
    }
  }
}





################################
# Early versions (do not revise)
################################

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
    
    colnames(res_out)<-c("Subgroup","n","n1","n0","m1","m0","RMST","HR (95% CI)","HR Adjusted")
    
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
    
    
    colnames(res)<-c("Subgroup","n","n1","n0","m1","m0","RMST","HR (95% CI)","HR Adjusted")
    
    return(list(res_out=res,rmst_fit=rmst_fit))
  }
}


plotband_survdiff<-function(res,xlab="Months"){
  plot.band(x=res$dpoints,mean=res$diff.dpoints,xlabel=xlab,
            ylabel=expression(delta(t)==hat(S)[1](t)-hat(S)[0](t)),
            band=TRUE,lower=res$sb.lower,upper=res$sb.upper,show.axes=TRUE,ltype="l")
  lines(res$dpoints,res$pw.upper,type="l",lty=2,lwd=0.50)
  lines(res$dpoints,res$pw.lower,type="l",lty=2,lwd=0.50)
  abline(h=0.0,lty=1,col="black",lwd=2)
  rug(res$dpoints)
}


# For plotting subgroups and the complement
# This is called withing forestsearch
# Quick and dirty plot
plot.subgroup<-function(tte.name,event.name,treat.name,wgt.name=NULL,sub1,sub1C,xloc1=NULL,xloc2=NULL,details=FALSE,show.logrank=FALSE,
                        ymin=0,cox.cex=0.7,prob.points=c(0,0.25,0.5,0.75,1.0),
                        cex_Yaxis=0.8,title1=NULL,title2=NULL,choose_ylim=TRUE, put.legend.arms='left',
                        exp.lab="Treat",con.lab="Control",legend.cex=0.70,risk.cex=0.70,yloc1=0.6,yloc2=0.6,subid=NULL,byrisk=2,fix.rows=TRUE,show.med=FALSE,
                        xlab="Months", ylab="Survival"){
  
  
  if(!is.null(subid) & (!is.null(title1) | !is.null(title2))) stop("subid OR titles need to be specified (not both)")
  
  # Set tpoints.add to span range of both subgroups
  aa<-max(sub1[,c(tte.name)])
  bb<-max(sub1C[,c(tte.name)])
  tpoints.add<-c(-1,max(aa,bb))
  
  if(fix.rows){
    old_par = par(no.readonly = TRUE)
    layout(matrix(1:2, ncol=2))
  }
  
  par(mar = old_par$mar + c(0, 0, 0, -2))
  
  
  df<-sub1
  
  km.fit <- KM.plot.2sample.weighted(df=df, tte.name=tte.name, event.name=event.name, treat.name=treat.name,
                                     tpoints.add=tpoints.add,
                                     risk.set=TRUE, by.risk=byrisk, risk.cex=risk.cex, censor.cex=risk.cex,
                                     risk_offset=0.15, risk_delta=0.075,
                                     Xlab=xlab,Ylab=ylab, cex_Yaxis=cex_Yaxis,
                                     show.ticks=TRUE,
                                     col.1="black", col.2="blue",
                                     arms = c(exp.lab,con.lab), 
                                     arm.cex=risk.cex, put.legend.arms=put.legend.arms,
                                     ltys=c(1,1),lwds=c(2,2),
                                     show.logrank=FALSE, lr.digits=2, put.legend.lr="top",
                                     qlabel="m =",
                                     show.med=TRUE, med.cex=risk.cex-0.15,
                                     show.cox=TRUE, cox.digits=3, cox.cex=risk.cex-0.10)
  
  
  
  
  if(!is.null(subid)) title(main=subid,cex.main=0.80)
  if(!is.null(title1)) title(main=title1,cex.main=0.80)
  
  par(mar = old_par$mar + c(0, -2, 0, 0))
  
  df<-sub1C
  
  km.fit <- KM.plot.2sample.weighted(df=df, tte.name=tte.name, event.name=event.name, treat.name=treat.name,
                                     tpoints.add=tpoints.add,
                                     risk.set=TRUE, by.risk=byrisk, risk.cex=risk.cex, censor.cex=risk.cex,
                                     risk_offset=0.15, risk_delta=0.075,
                                     Xlab=xlab,Ylab="", show.Y.axis=FALSE, cex_Yaxis=cex_Yaxis,
                                     show.ticks=TRUE,
                                     col.1="black", col.2="blue",
                                     show_arm_legend=FALSE, arms = c(exp.lab,con.lab), arm.cex=risk.cex,
                                     ltys=c(1,1),lwds=c(2,2),
                                     show.logrank=FALSE, lr.digits=2, put.legend.lr="top",
                                     qlabel="m =",
                                     show.med=TRUE, med.cex=risk.cex-0.15,
                                     show.cox=TRUE, cox.digits=3, cox.cex=risk.cex-0.10)
  
  if(!is.null(subid)) title(main="Complement",cex.main=0.80)
  if(!is.null(title2)) title(main=title2,cex.main=0.80)
  
  par(old_par)
  # Restore layout.
  layout(1:1)
  
}


# Forestploter format
SG_HRtable<-function(fs.est,fs_bc,treat.name,est.scale="hr",sg1_name,sg0_name,E.name="E",C.name="C"){

  ans <- list()
  
    resH <- fs_bc$H_estimates[,c("H2","H2_lower","H2_upper")]
    # Describing CIs
    Hstat<-c(unlist(resH[1,]))[c(1,2,3)]
    Hstat<-round(Hstat,3)
    # Return est,low,hi,se
    est <- Hstat[1]
    low <- Hstat[2]
    hi <- Hstat[3]
    se <- c((hi-est)/1.96)
    hr_H <- c(est,low,hi,se)
    #names(hr_H) <- c("est","low","hi","se")
    
    resH <- fs_bc$Hc_estimates[,c("H2","H2_lower","H2_upper")]
    # Describing CIs
    Hstat<-c(unlist(resH[1,]))[c(1,2,3)]
    Hstat<-round(Hstat,3)
    est <- Hstat[1]
    low <- Hstat[2]
    hi <- Hstat[3]
    se <- c((hi-est)/1.96)
    hr_Hc <- c(est,low,hi,se)
    #names(hr_Hc) <- c("est","low","hi","se")
    
    if(est.scale=="hr"){
    df1 <- subset(fs.est$df.est,treat.recommend==1)
    ntreat <- sum(df1[,treat.name])
    ncontrol <- sum(1-df1[,treat.name])
    
    hr_B <- c(ntreat,ncontrol,hr_Hc)
    #names(hr_B) <- c("E","C",names(hr_Hc))
    names(hr_B) <- c(E.name,C.name,"est","low","hi","se")
    
    aa <- as.data.frame(t(hr_B))
    aa$Subgroup <- c(sg1_name)
    
    df0 <- subset(fs.est$df.est,treat.recommend==0)
    ntreat <- sum(df0[,treat.name])
    ncontrol <- sum(1-df0[,treat.name])
    hr_notB <- c(ntreat,ncontrol,hr_H)
    #names(hr_notB) <- c("E","C",names(hr_H))
    names(hr_notB) <- c(E.name,C.name,"est","low","hi","se")
    
    bb <- as.data.frame(t(hr_notB))
    bb$Subgroup <- c(sg0_name)
    
    ans$out_B <- aa[,c("Subgroup",E.name,C.name,"est","low","hi","se")]
    ans$out_nonB <- bb[,c("Subgroup",E.name,C.name,"est","low","hi","se")]
    }
    
    if(est.scale=="1/hr"){
      df1 <- subset(fs.est$df.est,treat.recommend==0)
      ntreat <- sum(df1[,treat.name])
      ncontrol <- sum(1-df1[,treat.name])
      
      hr_B <- c(ntreat,ncontrol,hr_H)
      #names(hr_B) <- c("E","C",names(hr_H))
      names(hr_B) <- c(E.name,C.name,"est","low","hi","se")
      aa <- as.data.frame(t(hr_B))
      aa$Subgroup <- c(sg1_name)
      
      df0 <- subset(fs.est$df.est,treat.recommend==1)
      ntreat <- sum(df0[,treat.name])
      ncontrol <- sum(1-df0[,treat.name])
      hr_notB <- c(ntreat,ncontrol,hr_Hc)
      #names(hr_notB) <- c("E","C",names(hr_Hc))
      names(hr_notB) <- c(E.name,C.name,"est","low","hi","se")
      bb <- as.data.frame(t(hr_notB))
      bb$Subgroup <- c(sg0_name)
      
      ans$out_B <- aa[,c("Subgroup",E.name,C.name,"est","low","hi","se")]
      ans$out_nonB <- bb[,c("Subgroup",E.name,C.name,"est","low","hi","se")]
      }
        return(ans)
  }
    
 
SG_HRtable2<-function(dfa=NULL,df.OOB,fa_SG,ntreat_fa=NA,ncontrol_fa=NA,outcome.name,event.name,treat.name,sg_name,E.name="E",C.name="C"){

if(is.null(df.OOB) & is.null(fa_SG)){  
  # ITT analysis
  sf<- c("Surv(")
  sf <- paste0(sf,outcome.name,",")
  sf <- paste0(sf,event.name,") ~ ")
  trt <- paste0(treat.name)
  sf <- paste0(sf,trt,"")
  cox.formula <- as.formula(paste(sf))
  
  hr <- summary(coxph(cox.formula,data=dfa))$conf.int[c(1,3,4)]
  
# hr <- summary(coxph(Surv(os_time,os_event)~combo, data=dfa))$conf.int[c(1,3,4)]
  
  ntreat <- sum(dfa[,treat.name])
#  dtreat <- sum(dfa[,treat.name]*dfa[,event.name])
#  ntreat <- paste0(ntreat," (")
#  ntreat <- paste0(ntreat,dtreat)
#  ntreat <- paste0(ntreat,")")
  
  
  ncontrol <- sum(1-dfa[,treat.name])
#  dcontrol <- sum(c(1-dfa[,treat.name])*c(dfa[,event.name]))
#  ncontrol <- paste0(ncontrol," (")
#  ncontrol <- paste0(ncontrol,dcontrol)
#  ncontrol <- paste0(ncontrol,")")
  
#  est <- round(hr[1],2)
#  low <- round(hr[2],2)
#  hi <- round(hr[3],2)
#  se <- round(c((hi-est)/1.96),2)
 
    est <- hr[1]
    low <- hr[2]
    hi <- hr[3]
    se <- c((hi-est)/1.96)
  
   hr <- c(ntreat,ncontrol,est,low,hi,se)
  names(hr) <- c(E.name,C.name,"est","low","hi","se")
  aa <- as.data.frame(t(hr))
  aa$Subgroup <- c(sg_name)  
  aa <- aa[,c("Subgroup",E.name,C.name,"est","low","hi","se")]
return(aa)
}  
  
if(is.null(df.OOB) & !is.null(fa_SG)){  
  resSG <- fa_SG[,c("H2","H2_lower","H2_upper")]
  # Describing CIs
  Hstat<-c(unlist(resSG[1,]))[c(1,2,3)]
  Hstat<-round(Hstat,3)
  # Return est,low,hi,se
  est <- Hstat[1]
  low <- Hstat[2]
  hi <- Hstat[3]
  se <- c((hi-est)/1.96)
  hr_SG <- c(est,low,hi,se)
  hr_SG <- c(ntreat_fa,ncontrol_fa,hr_SG)
  names(hr_SG) <- c(E.name,C.name,"est","low","hi","se")
  aa <- as.data.frame(t(hr_SG))
  # Return sg_name if no OOB
  aa$Subgroup <- c(sg_name)
  aa <- aa[,c("Subgroup",E.name,C.name,"est","low","hi","se")]
  return(aa)
}
  
  if(!is.null(df.OOB) & !is.null(fa_SG)){  
    ans <- list()
    resSG <- fa_SG[,c("H2","H2_lower","H2_upper")]
    # Describing CIs
    Hstat<-c(unlist(resSG[1,]))[c(1,2,3)]
    Hstat<-round(Hstat,3)
    # Return est,low,hi,se
    est <- Hstat[1]
    low <- Hstat[2]
    hi <- Hstat[3]
    se <- c((hi-est)/1.96)
    hr_SG <- c(est,low,hi,se)
    hr_SG <- c(ntreat_fa,ncontrol_fa,hr_SG)
    names(hr_SG) <- c(E.name,C.name,"est","low","hi","se")
    aa <- as.data.frame(t(hr_SG))
    # Return full-analysis
    aa$Subgroup <- c("full-Analysis")
    aa <- aa[,c("Subgroup",E.name,C.name,"est","low","hi","se")]
    
    ans$fa <- aa

    sf<- c("Surv(")
    sf <- paste0(sf,outcome.name,",")
    sf <- paste0(sf,event.name,") ~ ")
    trt <- paste0(treat.name)
    sf <- paste0(sf,trt,"")
    cox.formula <- as.formula(paste(sf))
    
    hr <- summary(coxph(cox.formula,data=df.OOB))$conf.int[c(1,3,4)]
    
    #hr <- summary(coxph(Surv(os_time,os_event)~combo, data=df.OOB))$conf.int[c(1,3,4)]
    ntreat <- sum(df.OOB[,treat.name])
    ncontrol <- sum(1-df.OOB[,treat.name])
    est <- hr[1]
    low <- hr[2]
    hi <- hr[3]
    se <- c((hi-est)/1.96)
    hr <- c(ntreat,ncontrol,est,low,hi,se)
    names(hr) <- c(E.name,C.name,"est","low","hi","se")
    bb <- as.data.frame(t(hr))
    bb$Subgroup <- c("N-fold")  
    bb <- bb[,c("Subgroup",E.name,C.name,"est","low","hi","se")]
    
    ans$oob <- bb
  
return(ans)
  }
}  
  
 

# Post-hoc subgroups
# Return dataframe for forestplot
df_sgforest <- function(fs_OOB, df_sg, fs_bc, outcome.name, treat.name, sg1_name,sg0_name,E.name,C.name,est.scale="hr"){
ans <- list()  
if(est.scale=="hr"){  
# Benefitting subgroup (Hc)  
df.OOB <- subset(fs_OOB,treat.recommend==1)
ntreat_fa <- nrow(subset(df_sg,treat.recommend==1 & treat==1))
ncontrol_fa <- nrow(subset(df_sg,treat.recommend==1 & treat==0))

temp <- SG_HRtable2(df.OOB=df.OOB,fa_SG=fs_bc$Hc_estimates,ntreat_fa=ntreat_fa,ncontrol_fa=ncontrol_fa,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,
sg_name=sg1_name,E.name=E.name,C.name=C.name)

bb <- temp$fa
cc <- temp$oob
# Add header 
aa <- bb 
aa$Subgroup <- sg1_name
aa[,c(E.name,C.name)] <- ""
aa[,c("est","low","hi","se")] <- NA
ans$res_sg1 <- rbind(aa,bb,cc)
# Complement
df.OOB <- subset(fs_OOB,treat.recommend==0)
ntreat_fa <- nrow(subset(df_sg,treat.recommend==0 & treat==1))
ncontrol_fa <- nrow(subset(df_sg,treat.recommend==0 & treat==0))

temp <- SG_HRtable2(df.OOB=df.OOB,fa_SG=fs_bc$H_estimates,ntreat_fa=ntreat_fa,ncontrol_fa=ncontrol_fa,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,
                    sg_name=sg0_name,E.name=E.name,C.name=C.name)
bb <- temp$fa
cc <- temp$oob
# Add header 
aa <- bb 
aa$Subgroup <- sg0_name
aa[,c(E.name,C.name)] <- ""
aa[,c("est","low","hi","se")] <- NA
ans$res_sg0 <- rbind(aa,bb,cc)
}

if(est.scale=="1/hr"){  
    # Reverse the role of treatment recommend and H <---> Hc
    # Benefitting subgroup (H)  
    df.OOB <- subset(fs_OOB,treat.recommend==0)
    ntreat_fa <- nrow(subset(df_sg,treat.recommend==0 & treat==1))
    ncontrol_fa <- nrow(subset(df_sg,treat.recommend==0 & treat==0))
    # Here, treat.name is "treat2"
    temp <- SG_HRtable2(df.OOB=df.OOB,fa_SG=fs_bc$H_estimates,ntreat_fa=ntreat_fa,ncontrol_fa=ncontrol_fa,outcome.name=outcome.name,event.name=event.name,treat.name="treat2",
                        sg_name=sg1_name,E.name=E.name,C.name=C.name)
    
    bb <- temp$fa
    cc <- temp$oob
    # Add header 
    aa <- bb 
    aa$Subgroup <- sg1_name
    aa[,c(E.name,C.name)] <- ""
    aa[,c("est","low","hi","se")] <- NA
    ans$res_sg1 <- rbind(aa,bb,cc)
    
    # Complement
    df.OOB <- subset(fs_OOB,treat.recommend==1)
    ntreat_fa <- nrow(subset(df_sg,treat.recommend==1 & treat==1))
    ncontrol_fa <- nrow(subset(df_sg,treat.recommend==1 & treat==0))
    
    temp <- SG_HRtable2(df.OOB=df.OOB,fa_SG=fs_bc$Hc_estimates,ntreat_fa=ntreat_fa,ncontrol_fa=ncontrol_fa,outcome.name=outcome.name,event.name=event.name,treat.name="treat2",
                        sg_name=sg0_name,E.name=E.name,C.name=C.name)
    bb <- temp$fa
    cc <- temp$oob
    # Add header 
    aa <- bb 
    aa$Subgroup <- sg0_name
    aa[,c(E.name,C.name)] <- ""
    aa[,c("est","low","hi","se")] <- NA
    ans$res_sg0 <- rbind(aa,bb,cc)
  }
  
return(ans)
}
  



df_sgforest2 <- function(df_sg, fs_bc, outcome.name, treat.name, sg1_name,sg0_name,E.name,C.name,est.scale="hr"){
  #ans <- list()  
  if(est.scale=="hr"){  
    # Benefitting subgroup (Hc)  
    df.OOB <- NULL
    ntreat_fa <- nrow(subset(df_sg,treat.recommend==1 & treat==1))
    ncontrol_fa <- nrow(subset(df_sg,treat.recommend==1 & treat==0))
    
    temp <- SG_HRtable2(df.OOB=df.OOB,fa_SG=fs_bc$Hc_estimates,ntreat_fa=ntreat_fa,ncontrol_fa=ncontrol_fa,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,
                        sg_name=sg1_name,E.name=E.name,C.name=C.name)
    
   # bb <- temp
    # Add header 
   # aa <- bb 
   # aa$Subgroup <- sg1_name
   #  aa[,c(E.name,C.name)] <- ""
   # aa[,c("est","low","hi","se")] <- NA
    #ans$res_sg1 <- rbind(aa,bb)
    
    res_sg1 <- temp
    
    # Complement
    df.OOB <- NULL
    ntreat_fa <- nrow(subset(df_sg,treat.recommend==0 & treat==1))
    ncontrol_fa <- nrow(subset(df_sg,treat.recommend==0 & treat==0))
    
    temp <- SG_HRtable2(df.OOB=df.OOB,fa_SG=fs_bc$H_estimates,ntreat_fa=ntreat_fa,ncontrol_fa=ncontrol_fa,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,
                        sg_name=sg0_name,E.name=E.name,C.name=C.name)
    #bb <- temp$fa
    #cc <- temp$oob
    # Add header 
    #aa <- bb 
    #aa$Subgroup <- sg0_name
    #aa[,c(E.name,C.name)] <- ""
    #aa[,c("est","low","hi","se")] <- NA
    #ans$res_sg0 <- rbind(aa,bb,cc)
    
    res_sg0 <- temp
    
    ans <- rbind(res_sg1,res_sg0)
    
  }
  
  if(est.scale=="1/hr"){  
    # Reverse the role of treatment recommend and H <---> Hc
    # Benefitting subgroup (H)  

    df.OOB <- NULL
    
    ntreat_fa <- nrow(subset(df_sg,treat.recommend==0 & treat==1))
    ncontrol_fa <- nrow(subset(df_sg,treat.recommend==0 & treat==0))
    # Here, treat.name is "treat2"
    temp <- SG_HRtable2(df.OOB=df.OOB,fa_SG=fs_bc$H_estimates,ntreat_fa=ntreat_fa,ncontrol_fa=ncontrol_fa,outcome.name=outcome.name,event.name=event.name,treat.name="treat2",
                        sg_name=sg1_name,E.name=E.name,C.name=C.name)
    
    bb <- temp$fa
    cc <- temp$oob
    # Add header 
    aa <- bb 
    aa$Subgroup <- sg1_name
    aa[,c(E.name,C.name)] <- ""
    aa[,c("est","low","hi","se")] <- NA
    ans$res_sg1 <- rbind(aa,bb,cc)
    
    # Complement
    df.OOB <- NULL
    
    ntreat_fa <- nrow(subset(df_sg,treat.recommend==1 & treat==1))
    ncontrol_fa <- nrow(subset(df_sg,treat.recommend==1 & treat==0))
    
    temp <- SG_HRtable2(df.OOB=df.OOB,fa_SG=fs_bc$Hc_estimates,ntreat_fa=ntreat_fa,ncontrol_fa=ncontrol_fa,outcome.name=outcome.name,event.name=event.name,treat.name="treat2",
                        sg_name=sg0_name,E.name=E.name,C.name=C.name)
    bb <- temp$fa
    cc <- temp$oob
    # Add header 
    aa <- bb 
    aa$Subgroup <- sg0_name
    aa[,c(E.name,C.name)] <- ""
    aa[,c("est","low","hi","se")] <- NA
    ans$res_sg0 <- rbind(aa,bb,cc)
  }
  
  return(ans)
}


forestplot_baseline <- function(df,confounders.name,arrow_text=c("favors Treatment","Control"),E.name=c("Treat"),C.name=c("Control"),outcome.name,event.name,treat.name,footnote_text=NULL,
                                xlim = c(0.25, 2.0), ticks_at = c(0.65, 1.0, 1.3),title_text=NULL){
  
  title_text <- NULL
  
  # ITT 
  res_itt <- SG_HRtable2(dfa=df,df.OOB=NULL,fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,
                         sg_name="ITT",E.name=E.name,C.name=C.name)
  # Add separator header 
  zz <- res_itt  
  zz$Subgroup <- "                 "
  zz[,c(E.name,C.name)] <- ""
  zz[,c("est","low","hi","se")] <- NA
  
  
  dt <- rbind(res_itt,zz)
  
  for(sg_index in 1:length(confounders.name)){
    
    sg_name <- confounders.name[sg_index]
    
    cov <- df[,sg_name]
    cov_values <- sort(unique(cov),decreasing=TRUE)
    
    if(length(cov_values)>=4){
      cut <- c(quantile(cov,c(0.5)))
      
      thiscut <- paste0(sg_name,"<=",cut)  
      dfsg <- subset(df,eval(parse(text=thiscut)))
      
      resA <- SG_HRtable2(dfa=dfsg,df.OOB=NULL, sg_name=thiscut,
                          fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,E.name=E.name,C.name=C.name)
      
      
      thiscut <- paste0(sg_name,">",cut)  
      dfsg <- subset(df,eval(parse(text=thiscut)))
      
      resB <- SG_HRtable2(dfa=dfsg,df.OOB=NULL, sg_name=thiscut,
                          fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,E.name=E.name,C.name=C.name)
      dt <- rbind(dt,resA,resB,zz)
      rm("resA","resB")
    }
    
    if(length(cov_values)==3){
      cut <- c(quantile(cov,c(0.5)))
      
      thiscut <- paste0(sg_name,"==",cov_values[1])  
      dfsg <- subset(df,eval(parse(text=thiscut)))
      
      resA <- SG_HRtable2(dfa=dfsg,df.OOB=NULL, sg_name=thiscut,
                          fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,E.name=E.name,C.name=C.name)
      
      thiscut <- paste0(sg_name,"==",cov_values[2])  
      dfsg <- subset(df,eval(parse(text=thiscut)))
      resB <- SG_HRtable2(dfa=dfsg,df.OOB=NULL, sg_name=thiscut,
      fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,E.name=E.name,C.name=C.name)
      
      thiscut <- paste0(sg_name,"==",cov_values[3])  
      dfsg <- subset(df,eval(parse(text=thiscut)))
      resC <- SG_HRtable2(dfa=dfsg,df.OOB=NULL, sg_name=thiscut,
                          fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,E.name=E.name,C.name=C.name)
      
      
      dt <- rbind(dt,resA,resB,resC,zz)
      rm("resA","resB","resC")
    }
    
    
    
    if(length(cov_values)==2){
      
      thiscut <- paste0(sg_name,"==",cov_values[1])  
      dfsg <- subset(df,eval(parse(text=thiscut)))
      
      resA <- SG_HRtable2(dfa=dfsg,df.OOB=NULL, sg_name=thiscut,
                          fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,E.name=E.name,C.name=C.name)
      
      
      thiscut <- paste0(sg_name,"==",cov_values[2])  
      dfsg <- subset(df,eval(parse(text=thiscut)))
      
      resB <- SG_HRtable2(dfa=dfsg,df.OOB=NULL, sg_name=thiscut,
                          fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,E.name=E.name,C.name=C.name)
      dt <- rbind(dt,resA,resB,zz)
      rm("resA","resB")
    }
    
  }  
  
  sg_colors <- c("white")
  
  tm <- forest_theme(core=list(fg_params=list(hjust = 1, x = 0.9),
                               bg_params=list(fill = sg_colors)),
                     colhead=list(fg_params=list(hjust=0.5, x=0.5)),
                     footnote_gp = gpar(cex = 0.7, fontface = "italic", col = "blue"))
  
  dt$` ` <- paste(rep(" ", 20), collapse = " ")
  
  # Create a confidence interval column to display
  dt$`HR (95% CI)` <- ifelse(is.na(dt$se), "",
                             sprintf("%.2f (%.2f to %.2f)", dt$est, dt$low, dt$hi))
  
  
  p <- forest(dt[,c(1:3, 8:9)],
              title=title_text,
              est = as.numeric(dt$est),
              lower = as.numeric(dt$low), 
              upper = as.numeric(dt$hi),
              sizes = 0.4,
              ci_column = 4,
              ref_line = 1,
              arrow_lab = arrow_text,
              xlim = xlim,
              ticks_at = ticks_at,
              footnote = footnote_text, theme=tm)
  
  
  return(list(dt=dt,p=p))
}


# cv is % found
# B is sensitivity for benefitting
# Q for "quesionable"
sens_text <- function(ks_fold,est.scale){
  cv <- fs_kfold$find_summary["Any"]
  if(est.scale=="hr"){
    Q <- fs_kfold$sens_summary["sens_H"]
    B <- fs_kfold$sens_summary["sens_Hc"]
  }
  if(est.scale=="1/hr"){
    Q <- fs_kfold$sens_summary["sens_Hc"]
    B <- fs_kfold$sens_summary["sens_H"]
  }
  cv_text <- paste0("CV found = ",c(round(100*cv,0)),"%")
  aa <- paste0(round(100*B,0),"%,")
  bb <- paste0(round(100*Q,0),"%")
  sense_text <- paste("Agreement(+,-) = ", aa, bb, collapse=",")
  sg_text <- paste(cv_text,sense_text, sep=", ")
  return(sg_text)
}

