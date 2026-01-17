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


# Functions used in bootstrap operating characteristics
# Replicates coxph
getci_Cox<-function(df,est,se,target,alpha=0.025,digits=3){
  a<-df[est]
  b<-df[se]
  if(!is.numeric(target)){
    c<-df[target]
  }
  if(is.numeric(target)){
    c<-target
  }
  log_lb<-log(a)-qnorm(0.975)*(b/a)
  log_ub<-log(a)+qnorm(0.975)*(b/a)
  lb<-exp(log_lb)
  ub<-exp(log_ub)
  cov<-ifelse(c>=lb & c<=ub,1,0)
  return(list(lb=round(lb,3),ub=round(ub,3),cover=cov))
}
# Return avg, SD rounded to 3
SummaryStat<-function(df,name,sigdig=2,includeSD=FALSE,showSD=TRUE){
  if(is.data.table(df)){
    df<-as.data.frame(df)
  }
  # include "SD=" in output
  if(!includeSD){  
    temp<-na.omit(df[,c(name)])
    m<-round(mean(temp),sigdig)
    sig<-round(sqrt(var(temp)),sigdig)
    # Check if binary
    if(all(as.numeric(temp) %in% c(0,1))){
      sig<-round(sqrt(m*(1-m)/length(temp)),sigdig)
    }
    if(sigdig<=4 & m<=0.001) m<-"<0.001"
    if(sigdig<=4 & sig<=0.001) sig<-"<0.001"
    out<-paste0(m," (")
    out<-paste0(out,sig)
    out<-paste0(out,")")
  }
  if(includeSD){  
    temp<-na.omit(df[,c(name)])
    m<-round(mean(temp),sigdig)
    sig<-round(sqrt(var(temp)),sigdig)
    # Check if binary
    if(all(as.numeric(temp) %in% c(0,1))){
      sig<-round(sqrt(m*(1-m)/length(m)),sigdig)
    }
    if(sigdig<=4 & m<=0.001) m<-"<0.001"
    if(sigdig<=4 & sig<=0.001) sig<-"<0.001"
    out<-paste0(m," (SD=")
    out<-paste0(out,sig)
    out<-paste0(out,")")
  }
  if(!showSD) out<-c(m)  
  return(out)
}

DiffRate<-function(df1,df2,name,sigdig=2){
  if(is.data.table(df1)){
    df1<-as.data.frame(df1)
    df2<-as.data.frame(df2)
  }
  temp<-na.omit(df1[,c(name)])
  m1<-mean(temp)
  temp<-na.omit(df2[,c(name)])
  m2<-mean(temp)
  out<-round(100*(m2-m1),1)
  return(out)
}

pow_size<-function(df,minsize=0,sigdig=3){
  rej<-with(df,mean(rej12 & size.Hc>=minsize))
  return(round(rej,sigdig))
}

