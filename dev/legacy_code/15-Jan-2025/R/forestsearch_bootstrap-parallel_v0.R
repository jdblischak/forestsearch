
# Get operator (parallel or serial)
get_operator <- function() {
  is_par <- foreach::getDoParWorkers() > 1
  if (is_par) {
    res <- foreach::`%dopar%`
  } else {
    res <- foreach::`%do%`
  }
  res
}



# Forestsearch requires c("grf","policytree","data.table","randomForest","survival","plyr")
# Remove SPlit

#' Parallelize arbitrary user-supplied expression
#'
#' @param expr R expression.
#' @param n Number of iterations.
#' @param seed Initial random seed.
#' @param counter Character string. Name of the loop counter variable. Default is `"i"`.
#' @param export Exported variables to be used in `expr`.
bootPar <- function(expr, boots, seed, counter = "i", export) {
  `%op%` <- get_operator()
  expr <- substitute(expr)
  res <- foreach::foreach(
    i = seq_len(boots),
    .export = export,
    .packages = c("survival","grf","policytree","randomForest","data.table","plyr"),
    .combine = "rbind.fill",
    .errorhandling = "pass"
  ) %op% {
    assign(counter, value = i)
    set.seed(seed + i - 1)
    eval(expr)
  }
  res
}

# use rbind here
bootYstar <- function(expr, boots, seed, counter = "i", export) {
  `%op%` <- get_operator()
  expr <- substitute(expr)
  res <- foreach::foreach(
    i = seq_len(boots),
    .export = export,
    .packages = c("data.table"),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %op% {
    assign(counter, value = i)
    set.seed(seed + i - 1)
    eval(expr)
  }
  res
}


fun_arg_list_boot2<-c("dgm","maxFollow","sim_aftm1_gbsg","forestsearch","fs.estimates.out","vi.grf.min","dummy","subgroup.search",
                "subgroup.consistency","acm.disjctif","max.minutes","ztrail","one.zero","maxk","nmin.fs","d.min","hr.threshold",
                "NA.CHR.Weighted","R.Weighted","N.Weighted","hr.consistency","stop.threshold","fs.splits","m1.threshold","pconsistency.threshold",
                "confounders.name","outcome.name","id.name","treat.name",
                "cox.formula.sim","cox.formula.adj.sim",
                "event.name",
                "get.FG","FG.transform","LS.int.FG","vt.subg.harm.survival","vt.threshold","treat.threshold","maxdepth","n.min","ntree","vt.data2",
                "formatRCTDataset2","vt.estimates.out","label.analyses","grf.subg.harm.survival","grf.estimates.out","dmin.grf",
                "frac.tau","sg_focus","is.RCT",
                "get.FS","get.VT","get.GRF","grf.fs.subg","muC.adj","oc_analyses","N","split_method","quiet",
                "cox.formula.boot","forestsearch_bootstrap","oc_analyses_FSboot","fsBoot.estimates.out","boots","est.loghr",
                "get_Ystar","count.id","fsboot_forparallel","get_Cox_sg")


fun_arg_list_boot<-c("forestsearch","fs.estimates.out","vi.grf.min","dummy","subgroup.search",
                      "subgroup.consistency","acm.disjctif","max.minutes","ztrail","one.zero","maxk","nmin.fs","d.min","hr.threshold",
                      "NA.CHR.Weighted","R.Weighted","N.Weighted","hr.consistency","stop.threshold","fs.splits","m1.threshold","pconsistency.threshold",
                      "confounders.name","outcome.name","id.name","treat.name",
                      "event.name","sg_focus","pstop_futile","is.RCT",
                      "get.FG","FG.transform","LS.int.FG","vt.subg.harm.survival",
                      "formatRCTDataset2","grf.subg.harm.survival","grf.estimates.out","dmin.grf","n.min","frac.tau",
                      "grf.fs.subg","split_method","quiet",
                      "cox.formula.boot","forestsearch_bootstrap","oc_analyses_FSboot","fsBoot.estimates.out","boots",
                       "est.loghr","ci_est","get_dfRes","get_targetEst",
                      "get_Ystar","count.id","fsboot_forparallel","get_Cox_sg",
                     "df_boot_analysis","fs.est","H_obs","Hc_obs","nb_boots","est.scale",
                     "details","df.analysis","df_boot","id0","NN")


get_Ystar<-function(boot){
NN<-nrow(df_boot_analysis)
id0<-c(1:NN)
set.seed(8316951+boot*1000)
id_boot<-sample(id0,size=NN,replace=TRUE)
df_boot<-df_boot_analysis[id_boot,]
# re-name id.boot for bootstrap data
df_boot$id_boot<-c(1:nrow(df_boot))
Ystar<-c(unlist(lapply(df_boot_analysis$id,count.id,dfb=df_boot)))
return(c(Ystar))
}



fsboot_forparallel<-function(boot){

  NN<-nrow(df_boot_analysis)

  grp.consistency<-fs.est$grp.consistency  

  # Observed data evaluated at H:   hrH_obs
  # Observed data evaluated at H*   hrHstar_obs
  # Bootstrap data evaluated at H:  hrH_star 
  # Bootstrap data evaluated at H*: hrHstar_star
  
  id0<-c(1:NN)
  set.seed(8316951+boot*1000)
  id_boot<-sample(id0,size=NN,replace=TRUE)
  df_boot<-df_boot_analysis[id_boot,]
  # re-name id.boot for bootstrap data
  df_boot$id_boot<-c(1:nrow(df_boot))
 
  ############################################
  # Bootstrap data evaluated at H:  H_star 
  ############################################
  
  fitH_star<-get_Cox_sg(df_sg=subset(df_boot,treat.recommend==0),cox.formula=cox.formula.boot)
  H_star<-fitH_star$est_obs
  
  fitHc_star<-get_Cox_sg(df_sg=subset(df_boot,treat.recommend==1),cox.formula=cox.formula.boot)
  Hc_star<-fitHc_star$est_obs
  
  # Bias corrections
  H_biasadj_1<-H_biasadj_2<-NA
  Hc_biasadj_1<-Hc_biasadj_2<-NA
  tmins_search<-NA
  max_sg_est<-NA
  
  # Reset treat.recommend
  # df and df_boot have observed treat.recommend estimates
  # Bootstrap versions will create bootstrap version of treat.recommend
  # Can also just drop
  # Also remove initial confounders as these will be "re-assigned" for bootstrap
  # per lasso and/or grf (if employed)
  
  # Note: Need to drop initial confounders
  
  drop.vars<-c(fs.est$confounders.candidate,"treat.recommend")
  dfnew<-df_boot_analysis[ ,!(names(df_boot_analysis) %in% drop.vars)]
  dfnew_boot<-df_boot[ ,!(names(df_boot) %in% drop.vars)]

# Re-do grf (set to NULL) for bootstrap data
tempB<-try(forestsearch(df.analysis=dfnew_boot,
df.predict=dfnew,
is.RCT=fs.est$is.RCT,
Allconfounders.name=fs.est$Allconfounders.name,
max_n_confounders=fs.est$max_n_confounders,
vi.grf.min=fs.est$vi.grf.min,
conf_force=fs.est$conf_force,
use_lasso=fs.est$use_lasso, 
use_grf=fs.est$use_grf, 
use_grf_only=fs.est$use_grf_only,
dmin.grf=fs.est$dmin.grf, 
grf_res=NULL, grf_cuts=NULL, 
grf_depth=fs.est$grf_depth,
frac.tau=fs.est$frac.tau,
sg_focus=fs.est$sg_focus,
details=(boot <= 3),
outcome.name=fs.est$outcome.name,treat.name=fs.est$treat.name,
event.name=fs.est$event.name,id.name=fs.est$id.name,
n.min=fs.est$n.min,
hr.threshold=fs.est$hr.threshold,
hr.consistency=fs.est$hr.consistency,
stop.threshold=fs.est$stop.threshold,
pstop_futile=fs.est$pstop_futile,
fs.splits=fs.est$fs.splits,
d0.min=fs.est$d0.min,d1.min=fs.est$d1.min,
pconsistency.threshold=fs.est$pconsistency.threshold,
max.minutes=fs.est$max.minutes,
maxk=fs.est$maxk,plot.sg=FALSE),TRUE)

  # FS with error output NA's (set above)
  # FS done (w/o error) but NO sg found (replace NAs above with what is calculated)
  if(!inherits(tempB,"try-error") & is.null(tempB$sg.harm)){
    
    tmins_search<-tempB$find.grps$time_search
    
    if(!is.null(tempB$find.grps$max_sg_est)){
      max_sg_est<-tempB$find.grps$max_sg_est
    }
  }
  # FS done AND sg found
  if(!inherits(tempB,"try-error") & !is.null(tempB$sg.harm)){# Boot estimable & sg found 
    
    df_PredBoot<-tempB$df.predict # Observed data at bootstrap estimate of H
    
    dfboot_PredBoot<-tempB$df.est # Bootstrap data at H* (bootstrap estimate of H)
    
    max_sg_est<-c(as.numeric(tempB$find.grps$max_sg_est))
    tmins_search<-c(as.numeric(tempB$find.grps$time_search))
    
    rm("tempB")
    
    ##############################################
    # Observed data evaluated at H*   hrHstar_obs
    ##############################################
    fitHstar_obs<-get_Cox_sg(df_sg=subset(df_PredBoot,treat.recommend==0),cox.formula=cox.formula.boot)
    Hstar_obs<-fitHstar_obs$est_obs
    rm("fitHstar_obs")
    fitHstar_star<-get_Cox_sg(df_sg=subset(dfboot_PredBoot,treat.recommend==0),cox.formula=cox.formula.boot)
    Hstar_star<-fitHstar_star$est_obs
    rm("fitHstar_star")
    
    H_biasadj_1<-H_obs-(Hstar_star-Hstar_obs)
    H_biasadj_2<-2*H_obs-(H_star+Hstar_star-Hstar_obs)
    
    # Hc
    fitHcstar_obs<-get_Cox_sg(df_sg=subset(df_PredBoot,treat.recommend==1),cox.formula=cox.formula.boot)
    Hcstar_obs<-fitHcstar_obs$est_obs
    rm("fitHcstar_obs")
    
    fitHcstar_star<-get_Cox_sg(df_sg=subset(dfboot_PredBoot,treat.recommend==1),cox.formula=cox.formula.boot)
    Hcstar_star<-fitHcstar_star$est_obs
    rm("fitHcstar_star")
    
    Hc_biasadj_1 <- Hc_obs-(Hcstar_star-Hcstar_obs)
    Hc_biasadj_2 <- 2*Hc_obs-(Hc_star+Hcstar_star-Hcstar_obs)
    
  } # Bootstrap estimable
  
if(boot <= 3) cat("Bootstrap done, B=",c(boot),"\n")

  dfres<-data.table(H_biasadj_1,H_biasadj_2,
                    Hc_biasadj_1,Hc_biasadj_2,
                    tmins_search,max_sg_est)
  return(dfres)
}




# On original scale of est.loghr
# Does NOT convert 
ci_est<-function(x,sd,alpha=0.025,scale="hr",est.loghr=TRUE){
  c_alpha<-qnorm(1-alpha)
  c_low<-x-c_alpha*sd
  c_up<-x+c_alpha*sd
  est<-x
  if(scale=="hr" & est.loghr){
    # Only convert if log(hr) scale
    est<-exp(x)
    sd<-exp(x)*sd
    new_low<-exp(c_low)
    new_up<-exp(c_up)
  }
  if(scale=="1/hr" & est.loghr){
    est<-exp(-x)
    sd<-exp(-x)*sd
    new_low<-exp(-c_up)
    new_up<-exp(-c_low)
  }
  length<-new_up-new_low
  return(list(length=length,lower=new_low,upper=new_up,sd=sd,est=est))
}


get_dfRes<-function(Hobs,seHobs,H1_adj,H2_adj=NULL,ystar,cov_method="standard",cov_trim=0.0,est.scale="hr",est.loghr=TRUE){
# Un-adjusted
cest<-ci_est(x=Hobs,sd=seHobs,scale=est.scale,est.loghr=est.loghr)
H0_lower<-cest$lower
H0_upper<-cest$upper
sdH0<-cest$sd
H0<-cest$est
rm("cest")

est<-get_targetEst(x=H1_adj,ystar=ystar,cov_method=cov_method,cov_trim=cov_trim)
q<-est$target_est
se_new<-est$sehat_new
rm("est")

cest<-ci_est(x=q,sd=se_new,scale=est.scale,est.loghr=est.loghr)
H1_lower<-cest$lower
H1_upper<-cest$upper
sdH1<-cest$sd
H1<-cest$est
rm("cest","q","se_new")

outres<-data.table(H0,sdH0,H0_lower,H0_upper,
                   H1,sdH1,H1_lower,H1_upper)

if(!is.null(H2_adj)){
est<-get_targetEst(x=H2_adj,ystar=ystar,cov_method=cov_method,cov_trim=cov_trim)
q<-est$target_est
se_new<-est$sehat_new
rm("est")

cest<-ci_est(x=q,sd=se_new,scale=est.scale,est.loghr=est.loghr)
H2_lower<-cest$lower
H2_upper<-cest$upper
sdH2<-cest$sd
H2<-cest$est
rm("cest","q","se_new")

outres<-data.table(H0,sdH0,H0_lower,H0_upper,
                H1,sdH1,H1_lower,H1_upper,
                H2,sdH2,H2_lower,H2_upper)
}

return(outres)
}



fsBoot.parallel.out<-function(df,cox.formula.boot,FSboots,sim_number,est.scale="hr",est.loghr=TRUE){
  
  fit<-get_Cox_sg(df_sg=df,cox.formula=cox.formula.sim,est.loghr=est.loghr)
  H_itt<-fit$est_obs
  seH_itt<-fit$se_obs
  
  cest<-ci_est(x=H_itt,sd=seH_itt,scale=est.scale,est.loghr=est.loghr)
  H_itt<-cest$est
  seH_itt<-cest$sd
  H_itt_lower<-cest$lower
  H_itt_upper<-cest$upper
  rm("cest")
  rm("fit")
  
  fit<-get_Cox_sg(df_sg=df,cox.formula=cox.formula.adj.sim,est.loghr=est.loghr)
  H_adj_itt<-fit$est_obs["treat"]
  seH_adj_itt<-fit$se_obs["treat"]
  
  cest<-ci_est(x=H_adj_itt,sd=seH_adj_itt,scale=est.scale,est.loghr=est.loghr)
  H_adj_itt<-cest$est
  seH_adj_itt<-cest$sd
  H_adj_itt_lower<-cest$lower
  H_adj_itt_upper<-cest$upper
  rm("cest")
  rm("fit")
  
  # hatH(H true)
  fitH<-get_Cox_sg(df_sg=subset(df,flag.harm==1),cox.formula=cox.formula.boot,est.loghr=est.loghr)
  H_true<-fitH$est_obs
  seH_true<-fitH$se_obs
  
  if(est.scale=="hr" & est.loghr){
  H_true<-exp(H_true)
  seH_true<-H_true*seH_true
  }
  
  sizeH_true<-with(df,sum(flag.harm))
  propH_true<-with(df,mean(flag.harm))
  
  fitHc<-get_Cox_sg(df_sg=subset(df,flag.harm==0),cox.formula=cox.formula.boot,est.loghr=est.loghr)
  Hc_true<-fitHc$est_obs
  seHc_true<-fitHc$se_obs
  
  if(est.scale=="hr" & est.loghr){
  Hc_true<-exp(Hc_true)
  seHc_true<-Hc_true*seHc_true
  }  
  rm("fitH","fitHc")
  
  anyH_found<-any(names(df)=="treat.recommend")
  
  # Oracle versions
  hr_1_H<-with(subset(df,flag.harm==1),mean(h1.potential))
  hr_0_H<-with(subset(df,flag.harm==1),mean(h0.potential))
  H_oracle_causal<-hr_1_H/hr_0_H
  
  hr_1_Hc<-with(subset(df,flag.harm==0),mean(h1.potential))
  hr_0_Hc<-with(subset(df,flag.harm==0),mean(h0.potential))
  Hc_oracle_causal<-hr_1_Hc/hr_0_Hc
  
  if(est.loghr & est.scale=="loghr"){
    H_oracle_causal<-log(H_oracle_causal)
    Hc_oracle_causal<-log(Hc_oracle_causal)
  }

  
  if(!anyH_found){
    # Causal hat(H)
    hatH_causal<-NA
    # Causal hat(Hc)
    hatHc_causal<-NA
  }
  
  if(is.null(FSboots) | !anyH_found){
    # initiate as NA
    H_obs<-seH_obs<-H1.bc<-seH1.bc<-NA
    H2.bc<-seH2.bc<-NA
    Hc_obs<-seHc_obs<-Hc1.bc<-seHc1.bc<-NA
    Hc2.bc<-seHc2.bc<-NA
    H0_cover1<-H0_cover2<-H0_cover3<-H0_length<-H0_cover4<-NA
    H1_cover1<-H1_cover2<-H1_cover3<-H1_length<-H1_cover4<-NA
    H2_cover1<-H2_cover2<-H2_cover3<-H2_length<-H2_cover4<-NA
    Hc0_cover1<-Hc0_cover2<-Hc0_cover3<-Hc0_length<-Hc0_cover4<-NA
    Hc1_cover1<-Hc1_cover2<-Hc1_cover3<-Hc1_length<-Hc1_cover4<-NA
    Hc2_cover1<-Hc2_cover2<-Hc2_cover3<-Hc2_length<-Hc2_cover4<-NA
    
  }
  
  
  if(anyH_found){
    # Causal hat(H)
    hr_1_hatH<-with(subset(df,treat.recommend==0),mean(h1.potential))
    hr_0_hatH<-with(subset(df,treat.recommend==0),mean(h0.potential))
    hatH_causal<-hr_1_hatH/hr_0_hatH
    
    # Causal hat(Hc)
    hr_1_hatHc<-with(subset(df,treat.recommend==1),mean(h1.potential))
    hr_0_hatHc<-with(subset(df,treat.recommend==1),mean(h0.potential))
    hatHc_causal<-hr_1_hatHc/hr_0_hatHc
  
    # These are limiting ("marginal") true (oracle) Cox effects
    H_causal<-dgm$hr.H.true
    Hc_causal<-dgm$hr.Hc.true
  
    if(est.loghr & est.scale=="loghr"){
    hatH_causal<-log(hatH_causal)
    hatHc_causal<-log(hatHc_causal)
    H_causal<-log(dgm$hr.H.true)
    Hc_causal<-log(dgm$hr.Hc.true)
    }
      }
  
  
  if(!is.null(FSboots) & anyH_found){
    
    # Estimates are on est.scale (converted to HR if est.loghr=TRUE & est.scale="hr")
    
    dfH<-FSboots$H_estimates
   
    H_obs<-dfH$H0
    seH_obs<-dfH$sdH0
    
    H1.bc<-dfH$H1
    seH1.bc<-dfH$sdH1
  
    H2.bc<-dfH$H2
    seH2.bc<-dfH$sdH2
    
    # H_true is the oracle estimator based on observed data
    # hatH_causal is potential-outcome causal hr for estimated H-population
    # H_causal is the limiting marginal (super-population) of true subgroup (dgm$hr.H.true) "H fixed"
    
    # Un-adjusted CI
    cest<-ci_cover(lower=dfH$H0_lower,upper=dfH$H0_upper,target=c(H_true,hatH_causal,H_causal,H_oracle_causal))
    H0_cover1<-cest$cover[1]
    H0_cover2<-cest$cover[2]
    H0_cover3<-cest$cover[3]
    H0_cover4<-cest$cover[4]
    
    H0_length<-cest$LC
    rm("cest")
    
    cest<-ci_cover(lower=dfH$H1_lower,upper=dfH$H1_upper,target=c(H_true,hatH_causal,H_causal,H_oracle_causal))
    H1_cover1<-cest$cover[1]
    H1_cover2<-cest$cover[2]
    H1_cover3<-cest$cover[3]
    H1_cover4<-cest$cover[4]
    
    H1_length<-cest$LC
    rm("cest")
    
    cest<-ci_cover(lower=dfH$H2_lower,upper=dfH$H2_upper,target=c(H_true,hatH_causal,H_causal,H_oracle_causal))
    H2_cover1<-cest$cover[1]
    H2_cover2<-cest$cover[2]
    H2_cover3<-cest$cover[3]
    H2_cover4<-cest$cover[4]
    
    H2_length<-cest$LC
    rm("cest")
    rm("dfH")
    
    # Hc
    dfHc<-FSboots$Hc_estimates
    
    Hc_obs<-dfHc$H0
    seHc_obs<-dfHc$sdH0
    
    Hc1.bc<-dfHc$H1
    seHc1.bc<-dfHc$sdH1
    
    Hc2.bc<-dfHc$H2
    seHc2.bc<-dfHc$sdH2
    
    # Un-adjusted CI
    cest<-ci_cover(lower=dfHc$H0_lower,upper=dfHc$H0_upper,target=c(Hc_true,hatHc_causal,Hc_causal,Hc_oracle_causal))
    Hc0_cover1<-cest$cover[1]
    Hc0_cover2<-cest$cover[2]
    Hc0_cover3<-cest$cover[3]
    Hc0_cover4<-cest$cover[4]
    
    Hc0_length<-cest$LC
    rm("cest")
    
    cest<-ci_cover(lower=dfHc$H1_lower,upper=dfHc$H1_upper,target=c(Hc_true,hatHc_causal,Hc_causal,Hc_oracle_causal))
    Hc1_cover1<-cest$cover[1]
    Hc1_cover2<-cest$cover[2]
    Hc1_cover3<-cest$cover[3]
    Hc1_cover4<-cest$cover[4]
    
    Hc1_length<-cest$LC
    rm("cest")
    
    cest<-ci_cover(lower=dfHc$H2_lower,upper=dfHc$H2_upper,target=c(Hc_true,hatHc_causal,Hc_causal,Hc_oracle_causal))
    Hc2_cover1<-cest$cover[1]
    Hc2_cover2<-cest$cover[2]
    Hc2_cover3<-cest$cover[3]
    Hc2_cover4<-cest$cover[4]
    Hc2_length<-cest$LC
    rm("cest")
    rm("dfHc")
      }
  
  df.res<-data.table(sim_number,sizeH_true,propH_true,
                     H_true,Hc_true,
                     hatH_causal,hatHc_causal,
                     H_oracle_causal, Hc_oracle_causal,
                     H_obs,H1.bc,H2.bc,
                     Hc_obs,Hc1.bc,Hc2.bc,
                     seH_true,seHc_true,
                     seH_obs,seHc_obs,
                     seH1.bc,seH2.bc,
                     seHc1.bc,seHc2.bc,
                     H0_cover1,H0_cover2,H0_cover3,H0_cover4,
                     H0_length,
                     H1_cover1,H1_cover2,H1_cover3,H1_cover4,
                     H1_length,
                     H2_cover1,H2_cover2,H2_cover3,H2_cover4,
                     H2_length,
                     Hc0_cover1,Hc0_cover2,Hc0_cover3,Hc0_cover4,
                     Hc0_length,
                     Hc1_cover1,Hc1_cover2,Hc1_cover3,Hc1_cover4,
                     Hc1_length,
                     Hc2_cover1,Hc2_cover2,Hc2_cover3,Hc2_cover4,
                     Hc2_length,
                     H_itt,seH_itt,H_itt_lower,H_itt_upper,
                     H_adj_itt,seH_adj_itt,H_adj_itt_lower,H_adj_itt_upper)
  return(df.res)
}




