
# sqrt(N)*bias version
dfB_out_sqrtbias<-function(df,dgm,N){
  df<-na.omit(df)
  # H estimates
  # fixed
  target<-c(unlist(dgm["hr.H.true"]))
  # Oracle estimator
  temp<-getci_Cox(df=df,est="H_true",se="seH_true",target=unlist(dgm["hr.H.true"]))
  est_length<-c(unlist(temp$ub))-c(unlist(temp$lb))
  
  # Causal
    b_hr <- dgm$b_hr.true
  aa <- c(b_hr["treat"])
  # interaction term "gamma"
  bb <- c(b_hr["zh"])
  
  temp2<-getci_Cox(df=df,est="H_true",se="seH_true",target=exp(aa+bb))
  
  Hstats_a<-with(df,
                 c(mean(H_true),
                   sqrt(var(H_true)),mean(seH_true),
                   min(H_true),max(H_true),
                   0,
                   sqrt(N)*(mean(c(H_true-H_oracle_causal))),
                   sqrt(N)*((mean(H_true)-target)),
                   1,
                   mean(temp2$cover,na.rm=TRUE),
                   mean(temp$cover,na.rm=TRUE)))
  
  # Observed estimator
  Hstats_b<-with(df,c(mean(H_obs),
                      sqrt(var(H_obs)),mean(seH_obs),
                      min(H_obs), max(H_obs),
                      sqrt(N)*(mean(c(H_obs-H_true))),
                      sqrt(N)*(mean(c(H_obs-hatH_causal))),
                      sqrt(N)*((mean(H_obs)-target)),
                      mean(H0_cover1),
                      mean(H0_cover2),mean(H0_cover3)))
  
  # Bias-corrected
  Hstats_c<-with(df,c(mean(H2.bc),
                      sqrt(var(H2.bc)),mean(seH2.bc),
                      min(H2.bc),max(H2.bc),
                      sqrt(N)*(mean(c(H2.bc-H_true))),
                      sqrt(N)*(mean(c(H2.bc-hatH_causal))),
                      sqrt(N)*((mean(H2.bc)-target)),
                      mean(H2_cover1),
                      mean(H2_cover2),mean(H2_cover3)))
  
  
  # Hc estimates
  target<-c(unlist(dgm["hr.Hc.true"]))
  # Oracle estimator
  temp<-getci_Cox(df=df,est="Hc_true",se="seHc_true",target=unlist(dgm["hr.Hc.true"]))
  
  est_length<-c(unlist(temp$ub))-c(unlist(temp$lb))
  
  temp2<-getci_Cox(df=df,est="Hc_true",se="seHc_true",target=exp(aa))
  
  Hcstats_a<-with(df,
                  c(mean(Hc_true),
                    sqrt(var(Hc_true)),mean(seHc_true),
                    min(Hc_true),max(Hc_true),
                    0,
                    sqrt(N)*(mean(c(Hc_true-Hc_oracle_causal))),
                    sqrt(N)*((mean(Hc_true)-target)),
                    1,
                    mean(temp2$cover,na.rm=TRUE),
                    mean(temp$cover,na.rm=TRUE)))
  
  # Observed estimator
  Hcstats_b<-with(df,c(mean(Hc_obs),
                       sqrt(var(Hc_obs)),mean(seHc_obs),
                       min(Hc_obs),max(Hc_obs),
                       sqrt(N)*(mean(c(Hc_obs-Hc_true))),
                       sqrt(N)*(mean(c(Hc_obs-hatHc_causal))),
                       sqrt(N)*((mean(Hc_obs)-target)),
                       mean(Hc0_cover1),
                       mean(Hc0_cover2),mean(Hc0_cover3)))
  
  
  # Bias-corrected
  Hcstats_c<-with(df,c(mean(Hc2.bc),
                       sqrt(var(Hc2.bc)),mean(seHc2.bc),
                       min(Hc2.bc), max(Hc2.bc),
                       sqrt(N)*(mean(c(Hc2.bc-Hc_true))),
                       sqrt(N)*(mean(c(Hc2.bc-hatHc_causal))),
                       sqrt(N)*((mean(Hc2.bc)-target)),
                       mean(Hc2_cover1),
                       mean(Hc2_cover2),mean(Hc2_cover3)))
  
  
  Hstats<-rbind(Hstats_a,Hstats_b,Hstats_c)
  Hcstats<-rbind(Hcstats_a,Hcstats_b,Hcstats_c)
  
  Tdf<-as.data.frame(rbind(Hstats,Hcstats))
 
  
  rnH <- c("H-oracle","H-estimate","H-bcorrect")
  rnHc <- c("$H^{c}$-oracle","$H^{c}$-estimate","$H^{c}$-bcorrect")
  
  rownames(Tdf)<-c(rnH,rnHc)
  
  
  cnames <- c("Avg","SD","$\\widehat{SD}$",
              "min", "max",
              "$\\hat{b}^{oracle}$",            
              "$\\hat{b}^{\\ddagger}$",            
              "$b^{\\dagger}$",            
              "$\\hat{C}^{oracle}$",            
              "$\\hat{C}^{\\ddagger}$",            
              "$C^{\\dagger}$")
  
  
  colnames(Tdf) <- cnames
  
  return(Tdf)
}


# This is default
# H2 denotes the bootstrap modification of Harell's approach
# Quantities in df are calculated in fsBoot.parallel.out()
# Source file = forestsearch_bootstrap-parallel_v0.R

dfB_out<-function(df,dgm){
df<-na.omit(df)
# H estimates
# fixed
target<-c(unlist(dgm["hr.H.true"]))
# Oracle estimator
temp<-getci_Cox(df=df,est="H_true",se="seH_true",target=unlist(dgm["hr.H.true"]))
est_length<-c(unlist(temp$ub))-c(unlist(temp$lb))
# Causal
b_hr <- dgm$b_hr.true
aa <- c(b_hr["treat"])
# interaction term "gamma"
bb <- c(b_hr["zh"])

temp2<-getci_Cox(df=df,est="H_true",se="seH_true",target=exp(aa+bb))

Hstats_a<-with(df,
c(mean(H_true),
sqrt(var(H_true)),mean(seH_true),
min(H_true),max(H_true),
0,
100*(mean(c(H_true-H_oracle_causal)/H_oracle_causal)),
100*((mean(H_true)-target)/target),
mean(est_length),
1,
mean(temp2$cover,na.rm=TRUE),
mean(temp$cover,na.rm=TRUE)))

# Observed estimator
Hstats_b<-with(df,c(mean(H_obs),
sqrt(var(H_obs)),mean(seH_obs),
min(H_obs), max(H_obs),
100*(mean(c(H_obs-H_true)/H_true)),
100*(mean(c(H_obs-hatH_causal)/hatH_causal)),
100*((mean(H_obs)-target)/target),
mean(H0_length),
mean(H0_cover1),
mean(H0_cover2),mean(H0_cover3)))

# Bias-corrected
Hstats_c<-with(df,c(mean(H2.bc),
sqrt(var(H2.bc)),mean(seH2.bc),
min(H2.bc),max(H2.bc),
100*(mean(c(H2.bc-H_true)/H_true)),
100*(mean(c(H2.bc-hatH_causal)/hatH_causal)),
100*((mean(H2.bc)-target)/target),
mean(H2_length),
mean(H2_cover1),
mean(H2_cover2),mean(H2_cover3)))


# Hc estimates
target<-c(unlist(dgm["hr.Hc.true"]))
# Oracle estimator
temp<-getci_Cox(df=df,est="Hc_true",se="seHc_true",target=unlist(dgm["hr.Hc.true"]))

est_length<-c(unlist(temp$ub))-c(unlist(temp$lb))

temp2<-getci_Cox(df=df,est="Hc_true",se="seHc_true",target=exp(aa))

Hcstats_a<-with(df,
c(mean(Hc_true),
sqrt(var(Hc_true)),mean(seHc_true),
min(Hc_true),max(Hc_true),
0,
100*(mean(c(Hc_true-Hc_oracle_causal)/Hc_oracle_causal)),
100*((mean(Hc_true)-target)/target),
mean(est_length),
1,
mean(temp2$cover,na.rm=TRUE),
mean(temp$cover,na.rm=TRUE)))

# Observed estimator
Hcstats_b<-with(df,c(mean(Hc_obs),
sqrt(var(Hc_obs)),mean(seHc_obs),
min(Hc_obs),max(Hc_obs),
100*(mean(c(Hc_obs-Hc_true)/Hc_true)),
100*(mean(c(Hc_obs-hatHc_causal)/hatHc_causal)),
100*((mean(Hc_obs)-target)/target),
mean(Hc0_length),
mean(Hc0_cover1),
mean(Hc0_cover2),mean(Hc0_cover3)))


# Bias-corrected
Hcstats_c<-with(df,c(mean(Hc2.bc),
sqrt(var(Hc2.bc)),mean(seHc2.bc),
min(Hc2.bc), max(Hc2.bc),
100*(mean(c(Hc2.bc-Hc_true)/Hc_true)),
100*(mean(c(Hc2.bc-hatHc_causal)/hatHc_causal)),
100*((mean(Hc2.bc)-target)/target),
mean(Hc2_length),
mean(Hc2_cover1),
mean(Hc2_cover2),mean(Hc2_cover3)))


Hstats<-rbind(Hstats_a,Hstats_b,Hstats_c)
Hcstats<-rbind(Hcstats_a,Hcstats_b,Hcstats_c)

Tdf<-as.data.frame(rbind(Hstats,Hcstats))

rnH <- c("H-oracle","H-estimate","H-bcorrect")
rnHc <- c("$H^{c}$-oracle","$H^{c}$-estimate","$H^{c}$-bcorrect")

rownames(Tdf)<-c(rnH,rnHc)

cnames <- c("Avg","SD","$\\widehat{SD}$",
            "min", "max",
"$\\hat{b}^{oracle}$",            
"$\\hat{b}^{\\ddagger}$",            
"$b^{\\dagger}$", 
"Length",
"$\\hat{C}^{oracle}$",            
"$\\hat{C}^{\\ddagger}$",            
"$C^{\\dagger}$")

colnames(Tdf) <- cnames

return(Tdf)
}


df_Hplots_bycover<-function(df_res,dgm,return_plot=TRUE,qmax=c(1,1,1,1),alpha_sg=0.025){
  dfA<-as.data.frame(df_res)
  dfA<-na.omit(dfA)
  # Observed
  temp<-getci_Cox(df=na.omit(dfA),est="H_obs",se="seH_obs",target=unlist(dgm["hr.H.true"]),alpha=alpha_sg)
  Hobs_l<-unlist(temp$lb)
  Hobs_u<-unlist(temp$ub)
  Hobs_cover<-unlist(temp$cover)
  # Bias-corrected
  temp<-getci_Cox(df=na.omit(dfA),est="H2.bc",se="seH2.bc",target=unlist(dgm["hr.H.true"]),alpha=alpha_sg)
  Hbc_l<-unlist(temp$lb)
  Hbc_u<-unlist(temp$ub)
  Hbc_cover<-unlist(temp$cover)
  # True
  temp<-getci_Cox(df=na.omit(dfA),est="H_true",se="seH_true",target=unlist(dgm["hr.H.true"]),alpha=alpha_sg)
  H_l<-unlist(temp$lb)
  H_u<-unlist(temp$ub)
  H_cover<-unlist(temp$cover)
  nStar<-nrow(dfA)
  dH_plot <- data.frame(
    x=c(1:(nStar*3)), 
    value1=c(H_l,Hobs_l,Hbc_l), 
    value2=c(H_u,Hobs_u,Hbc_u),
    value3=c(H_cover,Hobs_cover,Hbc_cover),
    ests=c(rep('H_true',nStar),rep('H_un',nStar), rep('H_bc', nStar))
  )
  nullH.value<-dgm$hr.H.true
  # Hc
  # Observed
  temp<-getci_Cox(df=na.omit(dfA),est="Hc_obs",se="seHc_obs",target=unlist(dgm["hr.Hc.true"]),alpha=alpha_sg)
  Hobs_l<-unlist(temp$lb)
  Hobs_u<-unlist(temp$ub)
  Hobs_cover<-unlist(temp$cover)
  # Bias-corrected
  temp<-getci_Cox(df=na.omit(dfA),est="Hc2.bc",se="seHc2.bc",target=unlist(dgm["hr.Hc.true"]),alpha=alpha_sg)
  Hbc_l<-unlist(temp$lb)
  Hbc_u<-unlist(temp$ub)
  Hbc_cover<-unlist(temp$cover)
  # True
  temp<-getci_Cox(df=na.omit(dfA),est="Hc_true",se="seHc_true",target=unlist(dgm["hr.Hc.true"]),alpha=alpha_sg)
  H_l<-unlist(temp$lb)
  H_u<-unlist(temp$ub)
  H_cover<-unlist(temp$cover)
  
  dHc_plot <- data.frame(
    x=c(1:(nStar*3)), 
    value1=c(H_l,Hobs_l,Hbc_l), 
    value2=c(H_u,Hobs_u,Hbc_u),
    value3=c(H_cover,Hobs_cover,Hbc_cover),
    ests=c(rep('Hc_true',nStar),rep('Hc_un',nStar), rep('Hc_bc', nStar))
  )
  nullHc.value<-dgm$hr.Hc.true
  
  plot_H_Hc<-NULL
  if(return_plot){
    
    dfP<-subset(dH_plot,value3==0)
    
    # Modify ylim to 95% quantile
    maxy<-quantile(dfP$value2,qmax[1])
    miny<-min(dfP$value1)
    
    plot1<-ggplot(dfP) +
      geom_segment( aes(x=x, xend=x, y=value1, yend=value2), color="grey")+ 
      geom_point( aes(x=x, y=value1, color=factor(ests)), size=1 ) +
      geom_point( aes(x=x, y=value2, color=factor(ests)), size=1 ) +
      geom_hline(yintercept=nullH.value, linetype='dashed', color='black', size=0.5)+
      coord_flip(ylim=c(miny,maxy))+
      theme()+
      labs(title = "", colour="Est") +
      theme(legend.position = "top",
            legend.text = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    dfP<-subset(dH_plot,value3==1)
    
    maxy<-quantile(dfP$value2,qmax[2])
    miny<-min(dfP$value1)
    
    plot2<-ggplot(dfP) +
      geom_segment( aes(x=x, xend=x, y=value1, yend=value2), color="grey")+ 
      geom_point( aes(x=x, y=value1, color=factor(ests)), size=1 ) +
      geom_point( aes(x=x, y=value2, color=factor(ests)), size=1 ) +
      geom_hline(yintercept=nullH.value, linetype='dashed', color='black', size=0.5) +
      coord_flip(ylim=c(miny,maxy))+
      theme()+
      labs(title = "", colour="Est") +
      theme(legend.position = "top",
            legend.text = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    
    dfP<-subset(dHc_plot,value3==0)
    
    maxy<-quantile(dfP$value2,qmax[3])
    miny<-min(dfP$value1)
    
    plot3<-ggplot(dfP) +
      geom_segment( aes(x=x, xend=x, y=value1, yend=value2), color="grey")+ 
      geom_point( aes(x=x, y=value1, color=factor(ests)), size=1 ) +
      geom_point( aes(x=x, y=value2, color=factor(ests)), size=1 ) +
      geom_hline(yintercept=nullHc.value, linetype='dashed', color='black', size=0.5) +
      coord_flip(ylim=c(miny,maxy))+
      theme()+
      labs(title = "", colour="Est") +
      theme(legend.position = "top",
            legend.text = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    
    dfP<-subset(dHc_plot,value3==1)
    
    maxy<-quantile(dfP$value2,qmax[4])
    miny<-min(dfP$value1)
    
    plot4<-ggplot(dfP) +
      geom_segment( aes(x=x, xend=x, y=value1, yend=value2), color="grey")+ 
      geom_point( aes(x=x, y=value1, color=factor(ests)), size=1 ) +
      geom_point( aes(x=x, y=value2, color=factor(ests)), size=1 ) +
      geom_hline(yintercept=nullHc.value, linetype='dashed', color='black', size=0.5) +
      coord_flip(ylim=c(miny,maxy))+
      theme()+
      labs(title = "", colour="Est") +
      theme(legend.position = "top",
            legend.text = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
  plot_H_Hc<-plot_grid(plot1,plot2,plot3,plot4,labels=c("(a)","(b)","(c)","(d)"))
    
    
  }
  
  return(list(dH_plot=dH_plot,dHc_plot=dHc_plot,nullH.value=nullH.value,nullHc.value=nullHc.value,plot_H_Hc=plot_H_Hc,nStar=nStar))
}



df_Hplots<-function(df_res,dgm,return_plot=TRUE,qmax=c(1,1,1,1),alpha_sg=0.025){
  dfA<-as.data.frame(df_res)
  dfA<-na.omit(dfA)
  # Observed
  temp<-getci_Cox(df=na.omit(dfA),est="H_obs",se="seH_obs",target=unlist(dgm["hr.H.true"]),alpha=alpha_sg)
  Hobs_l<-unlist(temp$lb)
  Hobs_u<-unlist(temp$ub)
  Hobs_cover<-unlist(temp$cover)
  # Bias-corrected
  temp<-getci_Cox(df=na.omit(dfA),est="H2.bc",se="seH2.bc",target=unlist(dgm["hr.H.true"]),alpha=alpha_sg)
  Hbc_l<-unlist(temp$lb)
  Hbc_u<-unlist(temp$ub)
  Hbc_cover<-unlist(temp$cover)
  # True
  temp<-getci_Cox(df=na.omit(dfA),est="H_true",se="seH_true",target=unlist(dgm["hr.H.true"]),alpha=alpha_sg)
  H_l<-unlist(temp$lb)
  H_u<-unlist(temp$ub)
  H_cover<-unlist(temp$cover)
  nStar<-nrow(dfA)
  dH_plot <- data.frame(
    x=c(1:(nStar*3)), 
    value1=c(H_l,Hobs_l,Hbc_l), 
    value2=c(H_u,Hobs_u,Hbc_u),
    value3=c(H_cover,Hobs_cover,Hbc_cover),
    ests=c(rep('H_true',nStar),rep('H_un',nStar), rep('H_bc', nStar))
  )
  nullH.value<-dgm$hr.H.true
  # Hc
  # Observed
  temp<-getci_Cox(df=na.omit(dfA),est="Hc_obs",se="seHc_obs",target=unlist(dgm["hr.Hc.true"]),alpha=alpha_sg)
  Hobs_l<-unlist(temp$lb)
  Hobs_u<-unlist(temp$ub)
  Hobs_cover<-unlist(temp$cover)
  # Bias-corrected
  temp<-getci_Cox(df=na.omit(dfA),est="Hc2.bc",se="seHc2.bc",target=unlist(dgm["hr.Hc.true"]),alpha=alpha_sg)
  Hbc_l<-unlist(temp$lb)
  Hbc_u<-unlist(temp$ub)
  Hbc_cover<-unlist(temp$cover)
  # True
  temp<-getci_Cox(df=na.omit(dfA),est="Hc_true",se="seHc_true",target=unlist(dgm["hr.Hc.true"]),alpha=alpha_sg)
  H_l<-unlist(temp$lb)
  H_u<-unlist(temp$ub)
  H_cover<-unlist(temp$cover)
  
  dHc_plot <- data.frame(
    x=c(1:(nStar*3)), 
    value1=c(H_l,Hobs_l,Hbc_l), 
    value2=c(H_u,Hobs_u,Hbc_u),
    value3=c(H_cover,Hobs_cover,Hbc_cover),
    ests=c(rep('Hc_true',nStar),rep('Hc_un',nStar), rep('Hc_bc', nStar))
  )
  nullHc.value<-dgm$hr.Hc.true
  
  plot_H_Hc<-NULL
  if(return_plot){
    
    dfP <- dH_plot
    
    # Modify ylim to 95% quantile
    maxy<-quantile(dfP$value2,qmax[1])
    miny<-min(dfP$value1)
    
    plot1<-ggplot(dfP) +
      geom_segment( aes(x=x, xend=x, y=value1, yend=value2), color="grey")+ 
      geom_point( aes(x=x, y=value1, color=factor(ests)), size=1 ) +
      geom_point( aes(x=x, y=value2, color=factor(ests)), size=1 ) +
      geom_hline(yintercept=nullH.value, linetype='dashed', color='black', size=0.5)+
      coord_flip(ylim=c(miny,maxy))+
      theme()+
      labs(title = "", colour="Est") +
      theme(legend.position = "top",
            legend.text = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    dfP <- dHc_plot
    
    maxy<-quantile(dfP$value2,qmax[3])
    miny<-min(dfP$value1)
    
    plot2<-ggplot(dfP) +
      geom_segment( aes(x=x, xend=x, y=value1, yend=value2), color="grey")+ 
      geom_point( aes(x=x, y=value1, color=factor(ests)), size=1 ) +
      geom_point( aes(x=x, y=value2, color=factor(ests)), size=1 ) +
      geom_hline(yintercept=nullHc.value, linetype='dashed', color='black', size=0.5) +
      coord_flip(ylim=c(miny,maxy))+
      theme()+
      labs(title = "", colour="Est") +
      theme(legend.position = "top",
            legend.text = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    
plot_H_Hc<-plot_grid(plot1,plot2,labels=c("(a)","(b)"))
    
    
  }
  
  return(list(dH_plot=dH_plot,dHc_plot=dHc_plot,nullH.value=nullH.value,nullHc.value=nullHc.value,plot_H_Hc=plot_H_Hc,nStar=nStar))
}

#order_analysis<-c("H true","FS","GRF.60","VT(24)","VT(36)")
#plot_breaks<-c("H true","FS","GRF.60","VT(24)","VT(36)")
#col_values<-c("brown","blue","orange","lightblue","darkgrey")
#details=TRUE
#dig_Hc=2
#dig_H=2
#hrH_this=2.5

# Power plots
df_PowerPlots<-function(flist,order_analysis=c("H true","FS","GRF.60","VT(24)","VT(36)"),
                        plot_breaks=c("H true","FS","GRF.60","VT(24)","VT(36)"),col_values=c("brown","blue","orange","lightblue","darkgrey"),
                        order_analysis_b=c("H true","FS","GRF.60","ITT"),
                        plot_breaks_b=c("H true","FS","GRF.60","ITT"),
                        col_values_b=c("brown","blue","orange","beige"),
                        details=TRUE,dig_Hc=2,dig_H=2,hrH_this=2.5,Hlegend_type=1){
  res_all<-NULL
  for(ff in 1:length(flist)){
    load(paste0(resdir,flist[ff]))
    # Add "oracle" analysis
    if(details)  cat("dgm: H,H-complement,ITT=",c(dgm$hr.H.true,dgm$hr.Hc.true,dgm$hr.causal),"\n")  
    res_oracle<-subset(res,analysis=="FS")  
    if(details)  cat("# of FS simulations=",c(nrow(res_oracle)),"\n")  
    res_GRF<-subset(res,analysis=="GRF.60")  
    if(details)  cat("# of GRF.60 simulations=",c(nrow(res_GRF)),"\n")  
    # Set estimates (oracle) to truths
    res_oracle$hr.H.hat<-res_oracle$hr.H.true
    res_oracle$l.H.hat<-res_oracle$l.H.true
    res_oracle$u.H.hat<-res_oracle$u.H.true
    res_oracle$hr.Hc.hat<-res_oracle$hr.Hc.true  
    res_oracle$l.Hc.hat<-res_oracle$l.Hc.true  
    res_oracle$u.Hc.hat<-res_oracle$u.Hc.true  
    res_oracle$analysis<-c("H true")
    # Add ITT
    res_itt<-subset(res,analysis=="FS")
    res_itt$hr.Hc.hat<-res_itt$hr.itt
    res_itt$analysis<-c("ITT")
    res_new<-rbind(res,res_oracle,res_itt)
    res_new$hrH<-round(dgm$hr.H.true,3)
    res_new$hrHc<-round(dgm$hr.Hc.true,3)
    res_all<-rbind(res_all,res_new)
  }
  
  if(details){
    print(table(res_all$analysis))
    print(names(res_all))
  }
  
  res_all<-within(res_all,{
    rej1<-ifelse(!is.na(l.H.hat) & l.H.hat>1.0,1,0)
    rej2<-ifelse(!is.na(u.Hc.hat) & u.Hc.hat<1.0,1,0)
    rej12<-rej1*rej2
    anyH<-ifelse(!is.na(hr.H.hat),1,0)
    hrH<-ifelse(is.na(hrH),0.0,hrH)
    ppvHc<-npv
    sensHc<-specificity
    ppvH<-ppv
    sensH<-sensitivity
    ppv_Hc<-factor(cut2.factor(npv,cutpoints=c(0,0.9,0.95,1),labels=c("<=90%","90-95%",">95%")))
    sens_Hc<-factor(cut2.factor(specificity,cutpoints=c(0,0.9,0.95,1),labels=c("<=90%","90-95%",">95%")))
    ppv_H<-factor(cut2.factor(ppv,cutpoints=c(0,0.9,0.95,1),labels=c("<=90%","90-95%",">95%")))
    sens_H<-factor(cut2.factor(sensitivity,cutpoints=c(0,0.9,0.95,1),labels=c("<=90%","90-95%",">95%")))
  })
  
  # Restrict to FS, GRF.60, and VT(24)
  dfPlots<-subset(res_all,analysis %in% order_analysis)
  
  dfPlots<-within(dfPlots,
                  {Hc_rel<-100*(hr.Hc.hat-hrHc)/hrHc
                  H_rel<-100*(hr.H.hat-hrH)/hrH
                  Est<-factor(analysis,levels=order_analysis)})
  
  # Remove H true analysis as classification metrics not applicable
  dfP_ppvHc<-subset(dfPlots,analysis!="H true") %>% 
    mutate(Est=factor(analysis,levels=order_analysis)) %>%
    mutate(ppvHc=as.factor(ppv_Hc)) %>%  
    group_by(Est,ppvHc) %>%
    summarize(Hc_rbias=mean(Hc_rel))
  
  dfP_sensHc<-subset(dfPlots,analysis!="H true") %>% 
    mutate(Est=factor(analysis,levels=order_analysis)) %>%
    mutate(sensHc=as.factor(sens_Hc)) %>%  
    group_by(Est,sensHc) %>%
    summarize(Hc_rbias=mean(Hc_rel))
  
  
  plot_sensHc<-subset(dfPlots, analysis!="H true") %>% ggplot(aes(sens_Hc,Hc_rel,color=Est))+
    geom_boxplot()+
    scale_color_manual(breaks=plot_breaks,values=col_values)+
    geom_hline(yintercept=c(0.0,10), linetype='dashed', color=c('black','grey'), linewidth=c(0.5,0.25))+
    # geom_text(aes(label=round(anyH_rate,3)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1)+
    #coord_flip()+
    theme()+
    labs(title = "", colour="Est") +
    labs(x="sens(Hc)", y="Hc relative bias")+
    theme(legend.position = "top",
          legend.text = element_text(size = 4))
  
  plot_ppvHc<-subset(dfPlots, analysis!="H true") %>% ggplot(aes(ppv_Hc,Hc_rel,color=Est))+
    geom_boxplot()+
    scale_color_manual(breaks=plot_breaks,values=col_values)+
    geom_hline(yintercept=c(0.0,10), linetype='dashed', color=c('black','grey'), linewidth=c(0.5,0.25))+
    # geom_text(aes(label=round(anyH_rate,3)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1)+
    #coord_flip()+
    theme()+
    labs(title = "", colour="Est") +
    labs(x="ppv(Hc)", y="Hc relative bias")+
    theme(legend.position = "top",
          legend.text = element_text(size = 4))
  
  # For H, restrict to "estimable H's"
  # Otherwise, NA's
  plot_sensH<-subset(dfPlots, analysis!="H true" & !is.na(hr.H.hat)) %>% ggplot(aes(sens_H,H_rel,color=Est))+
    geom_boxplot()+
    scale_color_manual(breaks=plot_breaks,values=col_values)+
    geom_hline(yintercept=c(0.0,10), linetype='dashed', color=c('black','grey'), linewidth=c(0.5,0.25))+
    # geom_text(aes(label=round(anyH_rate,3)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1)+
    #coord_flip()+
    theme()+
    labs(title = "", colour="Est") +
    labs(x="sens(H)", y="H relative bias")+
    theme(legend.position = "top",
          legend.text = element_text(size = 4))
  
  plot_ppvH<-subset(dfPlots, analysis!="H true" & !is.na(hr.H.hat)) %>% ggplot(aes(ppv_H,H_rel,color=Est))+
    geom_boxplot()+
    scale_color_manual(breaks=plot_breaks,values=col_values)+
    geom_hline(yintercept=c(0.0,10), linetype='dashed', color=c('black','grey'), linewidth=c(0.5,0.25))+
    # geom_text(aes(label=round(anyH_rate,3)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1)+
    #coord_flip()+
    theme()+
    labs(title = "", colour="Est") +
    labs(x="ppv(H)", y="H relative bias")+
    theme(legend.position = "top",
          legend.text = element_text(size = 4))
  
  
  dfP_0<-dfPlots %>% 
    mutate(Est=factor(analysis,levels=order_analysis)) %>%
    mutate(hrH=as.factor(hrH)) %>%  
    group_by(Est,hrH) %>%
    summarize(anyH_rate=mean(anyH))
  
  dfP_1<-dfPlots %>% 
    mutate(Est=factor(analysis,levels=order_analysis)) %>%
    mutate(hrH=as.factor(hrH)) %>%  
    group_by(Est,hrH) %>%
    summarize(rej_rate=mean(rej1))
  
  # Conditional on finding H
  dfPlots_H<-subset(dfPlots,!is.na(hr.H.hat))
  
  dfP_2<-dfPlots_H %>% 
    mutate(Est=factor(analysis,levels=order_analysis)) %>%
    mutate(hrH=as.factor(hrH)) %>%  
    group_by(Est,hrH) %>%
    summarize(rej_rate=mean(rej2))
  
  dfP_12<-dfPlots %>% 
    mutate(Est=factor(analysis,levels=order_analysis)) %>%
    mutate(hrH=as.factor(hrH)) %>%  
    group_by(Est,hrH) %>%
    summarize(rej_rate=mean(rej12))
  
  plot0<-dfP_0 %>% ggplot(aes(hrH,anyH_rate,fill=Est))+
    geom_col(position="dodge")+
    scale_fill_manual(breaks=plot_breaks,values=col_values)+
    geom_hline(yintercept=0.025, linetype='dashed', color='black', linewidth=0.5)+
    # geom_text(aes(label=round(anyH_rate,3)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1)+
    coord_flip()+
    theme()+
    labs(title = "", colour="Est") +
    labs(x="Hazard ratio for H", y="Chance of identifying any H")+
    theme(legend.position = "top",
          legend.text = element_text(size = 6))
  
  plot1<-dfP_1 %>% ggplot(aes(hrH,rej_rate,fill=Est))+
    geom_col(position="dodge")+
    scale_fill_manual(breaks=plot_breaks,values=col_values)+
    geom_hline(yintercept=0.025, linetype='dashed', color='black', linewidth=0.5)+
    #geom_text(aes(label=round(rej_rate,2)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1,size=2.5)+
    coord_flip()+
    theme()+
    labs(title = "", colour="Est") +
    labs(x="Hazard ratio for H", y="Power for H > 1")+
    theme(legend.position = "top",
          legend.text = element_text(size = 6))
  
  # Remove x-axis here
  plot2<-dfP_2 %>% ggplot(aes(hrH,rej_rate,fill=Est))+
    geom_col(position="dodge")+
    scale_fill_manual(breaks=plot_breaks,values=col_values)+
    geom_hline(yintercept=0.025, linetype='dashed', color='black', linewidth=0.5)+
    #geom_text(aes(label=round(rej_rate,2)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1,size=2.5)+
    coord_flip()+
    theme()+
    labs(title = "", colour="Est") +
    labs(x="Hazard ratio for H", y="Conditional power for H-complement < 1")+
    theme(legend.position = "top",
          legend.text = element_text(size = 6),
          axis.title.y = element_blank())
  
  plot12<-dfP_12 %>% ggplot(aes(hrH,rej_rate,fill=Est))+
    geom_col(position="dodge")+
    scale_fill_manual(breaks=plot_breaks,values=col_values)+
    geom_hline(yintercept=0.025, linetype='dashed', color='black', linewidth=0.5)+
    #geom_text(aes(label=round(rej_rate,2)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1,size=2.5)+
    coord_flip()+
    theme()+
    labs(title = "", colour="Est") +
    labs(x="Hazard ratio for H", y="Power for {H > 1 & Hc < 1}")+
    theme(legend.position = "top",
          legend.text = element_text(size = 5),
          axis.title.y = element_blank())
  
  plot_all<-plot_grid(plot1,plot2,plot12,labels=c("(a)","(b)","(c)"), ncol=2, byrow=TRUE, align="hv")
  
  
  res_all_b<-within(res_all,{
    analysis2<-factor(analysis,levels=order_analysis_b)})
  
  if(!any(res_all$analysis=="FS-M")){
    df_FS<-subset(res_all_b, analysis=="FS" & !is.na(l.Hc.hat) & hrH==hrH_this)
    df_GRF<-subset(res_all_b, analysis=="GRF.60" & !is.na(l.Hc.hat) & hrH==hrH_this)
    df_true<-subset(res_all_b, analysis=="H true" & !is.na(l.Hc.hat) & hrH==hrH_this)
    df_itt<-subset(res_all_b, analysis=="ITT" & hrH==hrH_this)
    df_3<-rbind(df_FS,df_GRF,df_true)
    # Include ITT only for complement
    df_3b<-rbind(df_FS,df_GRF,df_true,df_itt)
  }
  
  if(any(res_all$analysis=="FS-M")){
    df_FS<-subset(res_all_b, analysis=="FS" & !is.na(l.Hc.hat) & hrH==hrH_this)
    df_FS_M<-subset(res_all_b, analysis=="FS-M" & !is.na(l.Hc.hat) & hrH==hrH_this)
    df_GRF<-subset(res_all_b, analysis=="GRF.60" & !is.na(l.Hc.hat) & hrH==hrH_this)
    df_true<-subset(res_all_b, analysis=="H true" & !is.na(l.Hc.hat) & hrH==hrH_this)
    df_itt<-subset(res_all_b, analysis=="ITT" & hrH==hrH_this)
    df_3<-rbind(df_FS,df_FS_M,df_GRF,df_true)
    # Include ITT only for complement
    df_3b<-rbind(df_FS,df_FS_M,df_GRF,df_true,df_itt)
  }
  
  #table(df_3$analysis)
  
  # Change analysis name to match above
  df_3$Est<-df_3$analysis2
  df_3b$Est<-df_3b$analysis2
  
  if(details) print(table(df_3b$analysis))
  
  # Relative bias
  plot3<-ggplot(df_3b,aes(Hc_rel,fill=Est))+
    geom_density(alpha=0.7)+
    geom_vline(xintercept=c(0.0,7.5), linetype='dashed', color='black', linewidth=c(0.15,0.35))+
    scale_fill_manual(breaks=plot_breaks,values=col_values)+
    theme()+
    labs(title = "", colour="Est") +
    labs(x="% Relative bias for H-complement", y="Density")+
    theme(legend.position = "top",
          legend.text = element_text(size = 6), axis.title.y = element_blank())
  
  # Remove the legend text and rely on plot12 for legend description
  if(Hlegend_type==1){
    # Hc Bias
    plot4<-ggplot(df_3b,aes(hr.Hc.hat,fill=Est))+
      geom_density(alpha=0.7)+
      geom_vline(xintercept=c(round(dgm$hr.Hc.true,dig_Hc)), linetype='dashed', color='black', linewidth=c(0.35))+
      scale_fill_manual(breaks=plot_breaks_b,values=col_values_b)+
      theme()+
      labs(title = "", colour="Est") +
      labs(x="HR estimates (Hc and ITT)", y="Density")+
      theme(legend.position = "top", legend.key.size=unit(0.5,'cm'), legend.title=element_blank(),
            legend.text = element_blank(), axis.title.y = element_blank())
    
    
    plot4b<-ggplot(df_3,aes(hr.H.hat,fill=Est))+
      geom_density(alpha=0.7)+
      geom_vline(xintercept=hrH_this, linetype='dashed', color='black', linewidth=c(0.35))+
      scale_fill_manual(breaks=plot_breaks,values=col_values)+
      theme()+
      labs(title = "", colour="Est") +
      labs(x="HR estimates (H)", y="Density") + 
      theme(legend.position = "top", legend.key.size=unit(0.5,'cm'), legend.title=element_blank(),
            legend.text = element_blank(), axis.title.y = element_blank())
  }
  
  if(Hlegend_type==2){
    # Hc Bias
    plot4<-ggplot(df_3b,aes(hr.Hc.hat,fill=Est))+
      geom_density(alpha=0.7)+
      geom_vline(xintercept=c(round(dgm$hr.Hc.true,dig_Hc)), linetype='dashed', color='black', linewidth=c(0.35))+
      scale_fill_manual(breaks=plot_breaks_b,values=col_values_b)+
      theme()+
      labs(title = "", colour="Est") +
      labs(x="HR estimates (Hc and ITT)", y="Density")+
      theme(legend.position = "top", legend.text = element_text(size=6), axis.title.y = element_blank())
    
    
    plot4b<-ggplot(df_3,aes(hr.H.hat,fill=Est))+
      geom_density(alpha=0.7)+
      geom_vline(xintercept=hrH_this, linetype='dashed', color='black', linewidth=c(0.35))+
      scale_fill_manual(breaks=plot_breaks,values=col_values)+
      theme()+
      labs(title = "", colour="Est") +
      labs(x="HR estimates (H)", y="Density") + 
      theme(legend.position = "top", legend.text = element_text(size=6), axis.title.y = element_blank())
  }
  
  
  
  
  plot_all2<-plot_grid(plot1,plot2,plot12,plot4,labels=c("(a)","(b)","(c)","(d)"), ncol=2, byrow=TRUE, align="hv")
  
  #if(details) plot(plot_all2)
  
  plot_all3<-plot_grid(plot4b,plot4,plot1,plot12,labels=c("(a)","(b)","(c)","(d)"), ncol=2, byrow=TRUE, align="hv")
  
  if(details){
    plot(plot_all3)
    #cat("Hc bias: OLS vs size Hc + ppv +npv + specificity +  sensitivity","\n")
    #ols_Hcbias<-lm(Hc_rel~ size.Hc+ppv+npv+specificity+sensitivity,data=dfPlots)
    #print(summary(ols_Hcbias))
    #cat("Hc bias: OLS vs size H + ppv +npv + specificity +  sensitivity","\n")
    #ols_Hcbias<-lm(Hc_rel~ size.H+ppv+npv+specificity+sensitivity,data=dfPlots)
    #print(summary(ols_Hcbias))
    #ols_Hcbias<-lm(Hc_rel~ specificity,data=dfPlots)
    #print(summary(ols_Hcbias))
    
    #ols_Hbias<-lm(H_rel~ size.Hc+ppv+npv+specificity+sensitivity,data=dfPlots)
    #cat("H bias: OLS vs size Hc + ppv +npv + specificity +  sensitivity","\n")
    #print(summary(ols_Hbias))
    #cat("H bias: OLS vs size H + ppv +npv + specificity +  sensitivity","\n")
    #ols_Hbias<-lm(H_rel~ size.H+ppv+npv+specificity+sensitivity,data=dfPlots)
    #print(summary(ols_Hbias))
    #ols_Hbias<-lm(H_rel~sensitivity,data=dfPlots)
    #print(summary(ols_Hbias))
  }
  
  
  return(list(dgm=dgm,res_all=res_all,dfPlots=dfPlots,plot0=plot0,plot_all2=plot_all2,
              plot4b=plot4b,plot4=plot4,plot_sensHc=plot_sensHc,plot_ppvHc=plot_ppvHc,
              plot_sensH=plot_sensH,plot_ppvH=plot_ppvH,
              df_3=df_3,df_3b=df_3b,plot_all3=plot_all3,plot1=plot1,plot2=plot2,plot12=plot12,
              order_analysis=order_analysis,order_analysis_b=order_analysis_b))
}
