count.id<-function(x,dfb){sum(dfb$id==x)}

calc_cov<-function(x,Est){mean(c((x-mean(x,na.rm=TRUE))*Est),na.rm=TRUE)}

calc_cov_Robust<-function(x,Est){
xmat<-cbind(x,Est)
xmat<-na.omit(xmat)
cov.nnve(xmat,emconv=0.1)$cov[1,2]
}

get_targetEst<-function(x,ystar,cov_method="standard",cov_trim=0.0){
 mx<-mean(x,na.rm=TRUE)
 # if(length(x) != nrow(ystar)) stop("Error in get_targetEst, x and ystar dimensions do not match")
    xc<-c(x-mx)
    N<-ncol(ystar)
    # Ystar_mat is (B x N) matrix
    # Across columns of Ystar 
    
  if(cov_method=="standard" | cov_method=="nocorrect")  cov_i<-apply(ystar,2,calc_cov,Est=xc)
    
  if(cov_method=="robust"){
     id.keep<-which(!is.na(x))
     xx<-x[id.keep]
     yystar<-ystar[id.keep,]
     cov_i<-apply(yystar,2,calc_cov_Robust,Est=xx)
  }  
    
    # Implement correction
    if(cov_method!="nocorrect"){
    varhat<-N*mean(cov_i^2)
    seH<-sqrt(varhat)
    # Denominator in cov_i is B.eval:
    B.eval<-sum(!is.na(xc))
    # correction
    nb_ratio<-N/(B.eval^2)
    termc<-nb_ratio*B.eval*mean(xc^2,na.rm=TRUE,trim=cov_trim)
    
    if(varhat<=termc){
    varhat_new<-varhat
    seH_new<-sqrt(varhat)
    }
  
    if(varhat>termc){
    varhat_new<-varhat-termc
    seH_new<-sqrt(varhat_new)
    }
  
    
      }
    if(cov_method=="nocorrect"){
    varhat<-N*mean(cov_i^2,trim=cov_trim)
    seH<-sqrt(varhat)
    termc<-0.0  
    varhat_new<-varhat
    seH_new<-seH
      }
    
    out<-(list(target_est=mx,sehat=seH,sehat_new=seH_new,term_correct=termc,varhat=varhat_new))
  return(out)
}

get_Cox_sg<-function(df_sg,cox.formula,est.loghr=TRUE){
names_tocheck<-all.vars(cox.formula)
check<-unlist(lapply(names_tocheck,grep,names(df_sg),value=TRUE))
check2<-match(names_tocheck,check)
if(sum(!is.na(check2))!=length(names_tocheck)) stop("df_sg dataset NOT contain cox.formula variables")
# May want to consider t-intervals on "\hat\theta" per Chang & Hall (2015)
# If so, could also improve with martingale-resampling 
fit<-summary(coxph(cox.formula,data=df_sg))$coefficients
# log(hr) parameters
if(est.loghr){
  bhat<-c(fit[,"coef"])  
  est_obs<-bhat
  se_obs<-c(fit[,"se(coef)"])
}
# Otherwise, hr 
if(!est.loghr){
  bhat<-c(fit[,"coef"])  
  est_obs<-exp(bhat)
  sebhat<-c(fit[,"se(coef)"])
  se_obs<-est_obs*sebhat
  }
return(list(est_obs=est_obs,se_obs=se_obs))
}


ci_cover<-function(lower,upper,target=0){
  if(length(target)==1){
    cover<-ifelse(lower<=target & upper>=target,1,0)
    LC<-upper-lower
    return(list(cover=cover,LC=LC))
  }
if(length(target)>1){
LC<-upper-lower
covers<-NULL
for(tt in 1:length(target)){
cover<-ifelse(lower<=target[tt] & upper>=target[tt],1,0)
covers<-c(covers,cover)
}
return(list(cover=covers,LC=LC))
  }
  }

# 2 Bootstrap targets
# outputs 2 (est,CIs)
getCIs<-function(Q1,Q2,ystar){
# Target 1  
est<-get_targetEst(x=Q1,ystar=ystar)
# If est.loghr=TRUE then on log(HR) scale
q1<-est$target_est
#se1<-est$sehat
se1_new<-est$sehat_new
se1<-est$sehat
rm("est")
# use SE new
cest<-ci_est(x=q1,sd=se1_new)
H1_lower<-cest$lower
H1_upper<-cest$upper
rm("cest")

# Target 2  
est<-get_targetEst(x=Q2,ystar=ystar)
# Call this H.bc
q2<-est$target_est
se2<-est$sehat
se2_new<-est$sehat_new
rm("est")
cest<-ci_est(x=q2,sd=se2_new)
H2_lower<-cest$lower
H2_upper<-cest$upper
rm("cest")
return(list(q1=q1,se1=se1,se1_new=se1_new,H1_lower=H1_lower,H1_upper=H1_upper,
            q2=q2,se2=se2,se2_new=se2_new,H2_lower=H2_lower,H2_upper=H2_upper))
}


# Used in simulations where flag.harm (true H) is known

forestsearch_bootstrap<-function(df,fs.est=NULL,boots=10,print.level=0,split_method="Random",
cox.formula.boot=as.formula(paste("Surv(y.sim,event.sim)~treat")),
cox2.formula.boot=as.formula(paste("Surv(Y,Event)~Treat")),
est.loghr=TRUE){

t.start<-proc.time()[3]
  
# Check names of df and df_boot
names_tocheck<-all.vars(cox.formula.boot)
check<-unlist(lapply(names_tocheck,grep,names(df),value=TRUE))
check2<-match(names_tocheck,check)
if(sum(!is.na(check2))!=length(names_tocheck)) stop("Dataset df does NOT contain cox.formula.boot variables")
# Create "id" variable if not included
if(!any(names(df)=="id")) df$id<-c(1:nrow(df))

NN<-nrow(df)

if(is.null(fs.est$grp.consistency$sg.harm)) stop("Forest search subgroup not found")
 
  grp.consistency<-fs.est$grp.consistency  
  sg.harm.fs<-grp.consistency$sg.harm
  
  # Observed data evaluated at H:   H_obs
  # Observed data evaluated at H*   Hstar_obs
  # Bootstrap data evaluated at H:  H_star 
  # Bootstrap data evaluated at H*: Hstar_star
  
  names_H<-c("H_obs","Hstar_obs","H_star","Hstar_star")
  names_Hc<-c("Hc_obs","Hcstar_obs","Hc_star","Hcstar_star")
  
  fs_H_boots<-as.data.frame(matrix(NA,nrow=boots,ncol=length(names_H)))
  fs_Hc_boots<-as.data.frame(matrix(NA,nrow=boots,ncol=length(names_Hc)))
  
  names(fs_H_boots)<-c(names_H)
  names(fs_Hc_boots)<-c(names_Hc)
  
  # H observed estimates
  # On est.loghr scale (TRUE = log(HR))
  fitH<-get_Cox_sg(df_sg=subset(df,treat.recommend==0),cox.formula=cox.formula.boot,est.loghr=est.loghr)
  H_obs<-fitH$est_obs
  seH_obs<-fitH$se_obs
  
  fs_H_boots[,"H_obs"]<-c(H_obs)
  
  # Hc observed estimates
  fitHc<-get_Cox_sg(df_sg=subset(df,treat.recommend==1),cox.formula=cox.formula.boot,est.loghr=est.loghr)
  Hc_obs<-fitHc$est_obs
  seHc_obs<-fitHc$se_obs
  fs_Hc_boots[,"Hc_obs"]<-c(Hc_obs)
  
  rm("fitH","fitHc")
  
  Ystar_mat<-matrix(NA,ncol=nrow(df),nrow=boots)
  
  # Observed data evaluated at H:   H_obs
  # Observed data evaluated at H*   Hstar_obs
  # Bootstrap data evaluated at H:  H_star 
  # Bootstrap data evaluated at H*: Hstar_star
  
  #vi_fs_covs<-data.frame(confounders=confounders.name,var.importance.boots=0)
  
  boots_eval<-0.0
  id0<-c(1:NN)
  
  for(boot in 1:boots){
    set.seed(8316951+boot*1000)
    
    id_boot<-sample(id0,size=NN,replace=TRUE)
    df_boot<-df[id_boot,]
    
    # re-name id.boot for bootstrap data
    df_boot$id_boot<-c(1:nrow(df_boot))
    
    Ystar<-c(unlist(lapply(df$id,count.id,dfb=df_boot)))
    
    Ystar_mat[boot,]<-c(Ystar)
  
    ############################################
    # Bootstrap data evaluated at H:  hrH_star 
    ############################################
    
    fitH_star<-get_Cox_sg(df_sg=subset(df_boot,treat.recommend==0),cox.formula=cox.formula.boot,est.loghr=est.loghr)
    H_star<-fitH_star$est_obs
    fs_H_boots[boot,"H_star"]<-c(H_star) 
    rm("fitH_star")
    
    fitHc_star<-get_Cox_sg(df_sg=subset(df_boot,treat.recommend==1),cox.formula=cox.formula.boot,est.loghr=est.loghr)
    Hc_star<-fitHc_star$est_obs
    fs_Hc_boots[boot,"Hc_star"]<-c(Hc_star)
    rm("fitHc_star")
    # Reset treat.recommend
    # df and df_boot have observed treat.recommend estimates
    # Bootstrap versions will create bootstrap estimate of treat.recommend
    # Drop treat.recommend from datasets
    drop.vars<-c("treat.recommend")
    dfnew<-df[ ,!(names(df) %in% drop.vars)]
    dfnew_boot<-df_boot[ ,!(names(df_boot) %in% drop.vars)]
    
    # id.name = id_boot?
    tempB<-try(forestsearch(df=dfnew_boot,confounders.name=confounders.name,
    df.predict=dfnew,
    split_method=split_method,
    details=FALSE,
    outcome.name=outcome.name,treat.name=treat.name,event.name=event.name,
    id.name=id.name,
    n.min=nmin.fs,hr.threshold=hr.threshold,hr.consistency=hr.consistency,fs.splits=fs.splits,
    stop.threshold=stop.threshold,d0.min=d.min,d1.min=d.min,
    pconsistency.threshold=pconsistency.threshold,
    max.minutes=max.minutes,
    sg_focus=sg_focus,
    maxk=maxk,plot.sg=FALSE,vi.grf.min=vi.grf.min),TRUE)
    
    if(!inherits(tempB,"try-error") & !is.null(tempB$sg.harm)){# Boot estimable & sg found 
      
      boots_eval<-boots_eval+1.0
      
      df_PredBoot<-tempB$df.pred # Observed data at bootstrap estimate of H
      dfboot_PredBoot<-tempB$df.est # Bootstrap data at H* (bootstrap estimate of H)
      
      names_tocheck<-all.vars(cox2.formula.boot)
      check<-unlist(lapply(names_tocheck,grep,names(dfboot_PredBoot),value=TRUE))
      check2<-match(names_tocheck,check)
      if(sum(!is.na(check2))!=length(names_tocheck)) stop("Bootstrap dataset (dfboot_PredBoot) df does NOT contain cox2.formula.boot variables")
      
      sg.harm.boot<-tempB$sg.harm
      
      rm("tempB")
    
      ##############################################
      # Observed data evaluated at H*   Hstar_obs
      ##############################################
      
      fitHstar_obs<-get_Cox_sg(df_sg=subset(df_PredBoot,treat.recommend==0),cox.formula=cox.formula.boot,est.loghr=est.loghr)
      Hstar_obs<-fitHstar_obs$est_obs
      fs_H_boots[boot,"Hstar_obs"]<-c(Hstar_obs)
      rm("fitHstar_obs")
      
      fitHcstar_obs<-get_Cox_sg(df_sg=subset(df_PredBoot,treat.recommend==1),cox.formula=cox.formula.boot,est.loghr=est.loghr)
      Hcstar_obs<-fitHcstar_obs$est_obs
      fs_Hc_boots[boot,"Hcstar_obs"]<-c(Hcstar_obs)
      rm("fitHcstar_obs")
      
      ###############################################
      # Bootstrap data evaluated at H*: Hstar_star
      ###############################################
     
      # Note, use cox2.formula.boot
      fitHstar_star<-get_Cox_sg(df_sg=subset(dfboot_PredBoot,treat.recommend==0),cox.formula=cox2.formula.boot,est.loghr=est.loghr)
      Hstar_star<-fitHstar_star$est_obs
      fs_H_boots[boot,"Hstar_star"]<-c(Hstar_star)
      rm("fitHstar_star")
      
      fitHcstar_star<-get_Cox_sg(df_sg=subset(dfboot_PredBoot,treat.recommend==1),cox.formula=cox2.formula.boot,est.loghr=est.loghr)
      Hcstar_star<-fitHcstar_star$est_obs
      fs_Hc_boots[boot,"Hcstar_star"]<-c(Hcstar_star)
      rm("fitHcstar_star")
  
      rm("df_boot","df2","df2_boot")
    
      } # Bootstrap estimable
    
  }  #Bootstrap loop
# H statistics 
  # H 
  # un-adjusted (Observed) CI
  # Call this C1 
  cest<-ci_est(x=H_obs,sd=seH_obs)
  H0_lower<-cest$lower
  H0_upper<-cest$upper
  rm("cest")
  
  # Observed data evaluated at H:   H_obs
  # Observed data evaluated at H*   Hstar_obs
  # Bootstrap data evaluated at H:  H_star 
  # Bootstrap data evaluated at H*: Hstar_star
  
# 4 CIs
# 1:biasadj: Hstar_star-Hstar_obs
# Hstar_star: boot data at boot H
# H_star_obs: observed data at boot H
# biasadj=H_obs-(Hstar_star-Hstar_obs)
# 2 (H_star_star-Hstar_obs)+(H_star-H_obs)
# biasadj=2*H_obs-(Hstar_star-Hstar_obs+H_star)
# 1 and 2 are with NO trimming
# 3 and 4 are (1 and 2) with trimming
# getCIs implements HR conversion if est.loghr
  
H_biasadj_1<-with(fs_H_boots,H_obs-(Hstar_star-Hstar_obs))
H_biasadj_2<-with(fs_H_boots,2*H_obs-(H_star+Hstar_star-Hstar_obs))

temp<-getCIs(Q1=c(H_biasadj_1),Q2=c(H_biasadj_2),ystar=Ystar_mat)

  H1.bc<-temp$q1
  seH1.bc<-temp$se1
  seH1.bc_new<-temp$se1_new
  H1_lower<-temp$H1_lower
  H1_upper<-temp$H1_upper
  
  H2.bc<-temp$q2
  seH2.bc<-temp$se2
  seH2.bc_new<-temp$se2_new
  H2_lower<-temp$H2_lower
  H2_upper<-temp$H2_upper
  
  rm("temp")
  
  dfH_out<-data.table(H_obs,seH_obs,H0_lower,H0_upper,
  H1.bc,seH1.bc,seH1.bc_new,H1_lower,H1_upper,
  H2.bc,seH2.bc,seH2.bc_new,H2_lower,H2_upper)
  
  # Hc statistics 
  
  cest<-ci_est(x=Hc_obs,sd=seHc_obs)
  Hc0_lower<-cest$lower
  Hc0_upper<-cest$upper
  rm("cest")

  Hc_biasadj_1<-with(fs_Hc_boots,Hc_obs-(Hcstar_star-Hcstar_obs))
  Hc_biasadj_2<-with(fs_Hc_boots,2*Hc_obs-(Hc_star+Hcstar_star-Hcstar_obs))
  
  temp<-getCIs(Q1=c(Hc_biasadj_1),Q2=c(Hc_biasadj_2),ystar=Ystar_mat)
  
  Hc1.bc<-temp$q1
  seHc1.bc<-temp$se1
  seHc1.bc_new<-temp$se1_new
  Hc1_lower<-temp$H1_lower
  Hc1_upper<-temp$H1_upper
  
  Hc2.bc<-temp$q2
  seHc2.bc<-temp$se2
  seHc2.bc_new<-temp$se2_new
  Hc2_lower<-temp$H2_lower
  Hc2_upper<-temp$H2_upper
  
  rm("temp")
  dfHc_out<-data.table(Hc_obs,seHc_obs,Hc0_lower,Hc0_upper,
                      Hc1.bc,seHc1.bc,seHc1.bc_new,Hc1_lower,Hc1_upper,
                      Hc2.bc,seHc2.bc,seHc2.bc_new,Hc2_lower,Hc2_upper)
  
t.now<-proc.time()[3]
t.min<-(t.now-t.start)/60

boot_out<-list(fs_H_boots=fs_H_boots,fs_Hc_boots=fs_Hc_boots,time_mins=t.min,
dfH_out=dfH_out,dfHc_out=dfHc_out)

return(boot_out)
}

get_dfPlot<-function(res,dgm){
hrH_true<-dgm$hr.H.true
hrHc_true<-dgm$hr.Hc.true
if(est.loghr & est.scale=="loghr"){
  hrH_true<-log(dgm$hr.H.true)
  hrHc_true<-log(dgm$hr.Hc.true)
}
res_new<-within(res,{
  b1H_1<-100*(H1.bc-H_true)/H_true
  b1H_2<-100*(H1.bc-hatH_causal)/hatH_causal
  b1H_3<-100*(H1.bc-hrH_true)/hrH_true
  
  b2H_1<-100*(H2.bc-H_true)/H_true
  b2H_2<-100*(H2.bc-hatH_causal)/hatH_causal
  b2H_3<-100*(H2.bc-hrH_true)/hrH_true
  
  b0H_1<-100*(H_obs-H_true)/H_true
  b0H_2<-100*(H_obs-hatH_causal)/hatH_causal
  b0H_3<-100*(H_obs-hrH_true)/hrH_true
}
)

df_bc<-NULL
# hr(H_true) target
res_new<-data.table(res_new)
df_res<-res_new[,c("b1H_1")]
df_res$est<-"BC_1"
df_res$target<-"H_true"
names(df_res)[1]<-"bias"
df_bc<-rbind(df_bc,df_res)

df_res<-res_new[,c("b2H_1")]
df_res$est<-"BC_2"
df_res$target<-"H_true"
names(df_res)[1]<-"bias"
df_bc<-rbind(df_bc,df_res)

df_res<-res_new[,c("b0H_1")]
df_res$est<-"Obs"
df_res$target<-"H_true"
names(df_res)[1]<-"bias"
df_bc<-rbind(df_bc,df_res)

# hatH_causal target
df_res<-res_new[,c("b1H_2")]
df_res$est<-"BC_1"
df_res$target<-"H_po"
names(df_res)[1]<-"bias"
df_bc<-rbind(df_bc,df_res)

df_res<-res_new[,c("b2H_2")]
df_res$est<-"BC_2"
df_res$target<-"H_po"
names(df_res)[1]<-"bias"
df_bc<-rbind(df_bc,df_res)

df_res<-res_new[,c("b0H_2")]
df_res$est<-"Obs"
df_res$target<-"H_po"
names(df_res)[1]<-"bias"
df_bc<-rbind(df_bc,df_res)

# H_causal fixed
df_res<-res_new[,c("b1H_3")]
df_res$est<-"BC_1"
df_res$target<-"H_fixed"
names(df_res)[1]<-"bias"
df_bc<-rbind(df_bc,df_res)

df_res<-res_new[,c("b2H_3")]
df_res$est<-"BC_2"
df_res$target<-"H_fixed"
names(df_res)[1]<-"bias"
df_bc<-rbind(df_bc,df_res)

df_res<-res_new[,c("b0H_3")]
df_res$est<-"Obs"
df_res$target<-"H_fixed"
names(df_res)[1]<-"bias"
df_bc<-rbind(df_bc,df_res)
return(df_bc)
}

var_summary<-function(res){
df_var<-NULL
# BC1 SD and Avg(est(sd))
aa<-with(res,sqrt(var(H1.bc,na.rm=TRUE)))
bb<-with(res,mean(seH1.bc,na.rm=TRUE))
cc<-with(res,mean(H1.bc,na.rm=TRUE))
dd<-with(res,mean(H1_cover1,na.rm=TRUE))
ee<-with(res,mean(H1_cover2,na.rm=TRUE))
ff<-with(res,mean(H1_cover3,na.rm=TRUE))
gg<-with(res,mean(H1_length,na.rm=TRUE))
df_var<-rbind(df_var,c(cc,aa,bb,dd,ee,ff,gg))

aa<-with(res,sqrt(var(H2.bc,na.rm=TRUE)))
bb<-with(res,mean(seH2.bc,na.rm=TRUE))
cc<-with(res,mean(H2.bc,na.rm=TRUE))
dd<-with(res,mean(H2_cover1,na.rm=TRUE))
ee<-with(res,mean(H2_cover2,na.rm=TRUE))
ff<-with(res,mean(H2_cover3,na.rm=TRUE))
gg<-with(res,mean(H2_length,na.rm=TRUE))
df_var<-rbind(df_var,c(cc,aa,bb,dd,ee,ff,gg))

aa<-with(res,sqrt(var(Hc1.bc,na.rm=TRUE)))
bb<-with(res,mean(seHc1.bc,na.rm=TRUE))
cc<-with(res,mean(Hc1.bc,na.rm=TRUE))
dd<-with(res,mean(Hc1_cover1,na.rm=TRUE))
ee<-with(res,mean(Hc1_cover2,na.rm=TRUE))
ff<-with(res,mean(Hc1_cover3,na.rm=TRUE))
gg<-with(res,mean(Hc1_length,na.rm=TRUE))
df_var<-rbind(df_var,c(cc,aa,bb,dd,ee,ff,gg))

aa<-with(res,sqrt(var(Hc2.bc,na.rm=TRUE)))
bb<-with(res,mean(seHc2.bc,na.rm=TRUE))
cc<-with(res,mean(Hc2.bc,na.rm=TRUE))
dd<-with(res,mean(Hc2_cover1,na.rm=TRUE))
ee<-with(res,mean(Hc2_cover2,na.rm=TRUE))
ff<-with(res,mean(Hc2_cover3,na.rm=TRUE))
gg<-with(res,mean(Hc2_length,na.rm=TRUE))
df_var<-rbind(df_var,c(cc,aa,bb,dd,ee,ff,gg))

colnames(df_var)<-c("Est","SD","Est(SD)","C1","C2","C3","L")
rownames(df_var)<-c("H1_bc","H2_bc","Hc1_bc","Hc2_bc")
return(round(df_var,digits=3))
}



