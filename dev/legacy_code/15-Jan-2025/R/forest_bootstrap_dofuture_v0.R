

forestsearch_bootstrap_dofuture<-function(fs.est,nb_boots,details=FALSE,show_three=FALSE){
  ## Bootstrap Cox formula
  sf<- c("Surv(")
  sf <- paste0(sf,fs.est$outcome.name,",")
  sf <- paste0(sf,fs.est$event.name,") ~ ")
  trt <- paste0(fs.est$treat.name)
  sf <- paste0(sf,trt,"")
  
  est.scale <- fs.est$est.scale
  
  cox.formula.boot <- as.formula(paste(sf))
  
  est.loghr <- TRUE
  
  df_boot_analysis <- fs.est$df.est
  
  fitH<-get_Cox_sg(df_sg=subset(df_boot_analysis,treat.recommend==0),cox.formula=cox.formula.boot)
  H_obs<-fitH$est_obs # log(hr) scale
  seH_obs<-fitH$se_obs
  # Hc observed estimates
  fitHc<-get_Cox_sg(df_sg=subset(df_boot_analysis,treat.recommend==1),cox.formula=cox.formula.boot)
  Hc_obs<-fitHc$est_obs
  seHc_obs<-fitHc$se_obs
  rm("fitH","fitHc")
  
  NN<-nrow(df_boot_analysis)
  id0<-c(1:NN)
  
  tB.start<-proc.time()[3]
  
  # message for backends
  # borrowed from simtrials (sim_fixed_n)
  if (!is(plan(), "sequential")) {
    # future backend
    message("Using ", nbrOfWorkers(), " cores with backend ", attr(plan("list")[[1]], "class")[2])
  } else if (foreach::getDoParWorkers() > 1) {
    message("Using ", foreach::getDoParWorkers(), " cores with backend ", foreach::getDoParName())
    message("Warning: ")
    message("doFuture may exhibit suboptimal performance when using a doParallel backend.")
  } else {
    message("Backend uses sequential processing.")
  }
  
  Ystar_mat <- foreach(
    boot = seq_len(nb_boots),
    .options.future=list(seed=TRUE),
    .combine="rbind", 
    .errorhandling="pass"
  ) %dofuture% {
    set.seed(8316951+boot*100)
    id_boot<-sample(id0,size=NN,replace=TRUE)
    df_boot<-df_boot_analysis[id_boot,]
    # re-name id.boot for bootstrap data
    df_boot$id_boot<-c(1:nrow(df_boot))
    ystar<-c(unlist(lapply(df_boot_analysis$id,count.id,dfb=df_boot)))
    return(c(ystar))
  }
  if(details) cat("Done with Ystar_mat","\n")
  
  results <- foreach(
    boot = seq_len(nb_boots),
    .options.future=list(seed=TRUE),
    .combine="rbind", 
    .errorhandling="pass"
  ) %dofuture% {
    show3 <- FALSE
    if(show_three) show3 <- (boot <= 3)
    # These bootstraps need to be aligned with ystar 
    # which is run separately
    # Using same seeds (as in ystar) to connect
    set.seed(8316951+boot*100)
    id_boot<-sample(id0,size=NN,replace=TRUE)
    df_boot<-df_boot_analysis[id_boot,]
    # re-name id.boot for bootstrap data
    df_boot$id_boot<-c(1:nrow(df_boot))
    ############################################
    # Bootstrap data evaluated at H:  H_star 
    ############################################
    fitH_star<-get_Cox_sg(df_sg=subset(df_boot,treat.recommend==0),cox.formula=cox.formula.boot, est.loghr=TRUE)
    H_star<-fitH_star$est_obs
    fitHc_star<-get_Cox_sg(df_sg=subset(df_boot,treat.recommend==1),cox.formula=cox.formula.boot, est.loghr=TRUE)
    Hc_star<-fitHc_star$est_obs
    
    # Bias corrections
    # Return NA if not estimable
    H_biasadj_1<-H_biasadj_2<-NA
    Hc_biasadj_1<-Hc_biasadj_2<-NA
    tmins_search<-NA
    max_sg_est<-NA
    prop_maxk <- NA
    L <- NA
    max_count <- NA
    
    # Reset treat.recommend
    # df and df_boot have observed treat.recommend estimates
    # Bootstrap versions will create bootstrap version of treat.recommend
    # Can also just drop
    # Also remove initial confounders as these will be re-assessed for bootstrap
    # per lasso and/or grf (if employed)
    
    # Note: Need to drop initial confounders
    drop.vars<-c(fs.est$confounders.candidate,"treat.recommend")
    dfnew<-df_boot_analysis[ ,!(names(df_boot_analysis) %in% drop.vars)]
    dfnew_boot<-df_boot[ ,!(names(df_boot) %in% drop.vars)]
    
    # Re-do grf (set to NULL) for bootstrap data
    # Note: est.scale is NOT imported here --> 
    # The re-scaling is done "at the end" below
    
    tempB<-try(forestsearch(df.analysis=dfnew_boot,
    df.predict=dfnew,
    is.RCT=fs.est$is.RCT,
    defaultcut_names=fs.est$defaultcut_names, cut_type=fs.est$cut_type,
    confounders.name=fs.est$confounders.name,
    max_n_confounders=fs.est$max_n_confounders,
    vi.grf.min=fs.est$vi.grf.min,
    conf_force=fs.est$conf_force,
    exclude_cuts=fs.est$exclude_cuts,
    replace_med_grf=fs.est$replace_med_grf,
    use_lasso=fs.est$use_lasso,use_grf=fs.est$use_grf,use_grf_only=fs.est$use_grf_only,
    dmin.grf=fs.est$dmin.grf,grf_res=NULL, grf_cuts=NULL,grf_depth=fs.est$grf_depth,
    frac.tau=fs.est$frac.tau,sg_focus=fs.est$sg_focus,
    details=show3,outcome.name=fs.est$outcome.name,treat.name=fs.est$treat.name,
    event.name=fs.est$event.name,id.name=fs.est$id.name,
    n.min=fs.est$n.min,hr.threshold=fs.est$hr.threshold,hr.consistency=fs.est$hr.consistency,
    stop.threshold=fs.est$stop.threshold,pstop_futile=fs.est$pstop_futile,fs.splits=fs.est$fs.splits,
    d0.min=fs.est$d0.min,d1.min=fs.est$d1.min,pconsistency.threshold=fs.est$pconsistency.threshold,
    max.minutes=fs.est$max.minutes,
    showten_subgroups=FALSE,
    maxk=fs.est$maxk,plot.sg=FALSE),TRUE)
    
    
    # FS with error output NA's (set above)
    # FS done (w/o error) but NO sg found (replace NAs above with what is calculated)
    if(!inherits(tempB,"try-error")){
    if(is.null(tempB$sg.harm)){
      tmins_search<-tempB$find.grps$time_search
      prop_maxk <- tempB$prop_maxk
      max_count <- tempB$find.grps$max_count
      L <- tempB$find.grps$L
      
      if(!is.null(tempB$find.grps$max_sg_est)){
        max_sg_est<-tempB$find.grps$max_sg_est
      }
    }
    # FS done AND sg found
    if(!is.null(tempB$sg.harm)){# Boot estimable & sg found 
      df_PredBoot<-tempB$df.predict # Observed data at bootstrap estimate of H
      dfboot_PredBoot<-tempB$df.est # Bootstrap data at H* (bootstrap estimate of H)
      max_sg_est<-c(as.numeric(tempB$find.grps$max_sg_est))
      tmins_search<-c(as.numeric(tempB$find.grps$time_search))
      prop_maxk <- c(as.numeric(tempB$prop_maxk))
      max_count <- tempB$find.grps$max_count
      L <- tempB$find.grps$L
      
      rm("tempB")
      ##############################################
      # Observed data evaluated at H*   hrHstar_obs
      ##############################################
      fitHstar_obs<-get_Cox_sg(df_sg=subset(df_PredBoot,treat.recommend==0),cox.formula=cox.formula.boot,est.loghr=TRUE)
      Hstar_obs<-fitHstar_obs$est_obs
      rm("fitHstar_obs")
      fitHstar_star<-get_Cox_sg(df_sg=subset(dfboot_PredBoot,treat.recommend==0),cox.formula=cox.formula.boot,est.loghr=TRUE)
      Hstar_star<-fitHstar_star$est_obs
      rm("fitHstar_star")
      
      H_biasadj_1<-H_obs-(Hstar_star-Hstar_obs)
      H_biasadj_2<-2*H_obs-(H_star+Hstar_star-Hstar_obs)
      
      # Hc
      fitHcstar_obs<-get_Cox_sg(df_sg=subset(df_PredBoot,treat.recommend==1),cox.formula=cox.formula.boot,est.loghr=TRUE)
      Hcstar_obs<-fitHcstar_obs$est_obs
      rm("fitHcstar_obs")
      
      fitHcstar_star<-get_Cox_sg(df_sg=subset(dfboot_PredBoot,treat.recommend==1),cox.formula=cox.formula.boot,est.loghr=TRUE)
      Hcstar_star<-fitHcstar_star$est_obs
      rm("fitHcstar_star")
      
      Hc_biasadj_1 <- Hc_obs-(Hcstar_star-Hcstar_obs)
      Hc_biasadj_2 <- 2*Hc_obs-(Hc_star+Hcstar_star-Hcstar_obs)
      
    } # Bootstrap estimable
    }
    
    if(show3) cat("Bootstrap done, B=",c(boot),"\n")
    
    dfres<-data.table(H_biasadj_1,H_biasadj_2,
                      Hc_biasadj_1,Hc_biasadj_2,
                      tmins_search,max_sg_est,prop_maxk,L,max_count)
    return(dfres)
  }
  
  tB.now<-proc.time()[3]
  tB.min<-(tB.now-tB.start)/60
  if(details){
    cat("Minutes for Boots",c(nb_boots,tB.min),"\n")
    cat("Projection per 1000",c(tB.min*(1000/nb_boots)),"\n")
    cat("Propn bootstrap subgroups found =",c(sum(!is.na(results$H_biasadj_2))/nb_boots),"\n")
    # How many timmed out 
    cat("Number timmed out=",c(sum(is.na(results$H_biasadj_2) & results$tmins_search>fs.est$max.minutes)),"\n")
    cat("Avg proportion of maxk searched",mean(results$prop_maxk),"\n")
  }
  
  # H and Hc are rescaled according to est.scale (ci_est)
  
  H_estimates<-try(get_dfRes(Hobs=H_obs,seHobs=seH_obs,H1_adj=results$H_biasadj_1,H2_adj=results$H_biasadj_2,
                             ystar=Ystar_mat,cov_method="standard",cov_trim=0.0,est.scale=est.scale,est.loghr=TRUE),TRUE)
  Hc_estimates<-try(get_dfRes(Hobs=Hc_obs,seHobs=seHc_obs,H1_adj=results$Hc_biasadj_1,H2_adj=results$Hc_biasadj_2,
                          ystar=Ystar_mat,cov_method="standard",cov_trim=0.0, est.scale=est.scale,est.loghr=TRUE),TRUE)
  
  if(inherits(H_estimates,"try-error") | inherits(Hc_estimates,"try-error")){
  out <- list(results=results,SG_CIs=NULL,FSsg_tab=NULL,Ystar_mat=Ystar_mat,H_estimates=NULL,Hc_estimates=NULL)
  }
  
  if(!inherits(H_estimates,"try-error") & !inherits(Hc_estimates,"try-error")){
  resH <- H_estimates[,c("H0","H0_lower","H0_upper")]
  # Describing CIs
  Hstat<-c(unlist(resH[1,]))[c(1,2,3)]
  Hstat<-round(Hstat,2)
  #a<-paste0(Hstat[1]," (95% CI=")
  a<-paste0(Hstat[1]," (")
  a<-paste0(a,Hstat[2])
  a<-paste0(a,",")
  a<-paste0(a,Hstat[3])
  H_res1<-paste0(a,")")
  ## Note: Ignore "H1 versions" 
  ## These are original Harrell versions but 
  ## does not perform well, here
  resH <- H_estimates[,c("H2","H2_lower","H2_upper")]
  # Describing CIs
  Hstat<-c(unlist(resH[1,]))[c(1,2,3)]
  Hstat<-round(Hstat,2)
  a<-paste0(Hstat[1]," (")
  a<-paste0(a,Hstat[2])
  a<-paste0(a,",")
  a<-paste0(a,Hstat[3])
  H_res2<-paste0(a,")")
  # Hc
  resH <- Hc_estimates[,c("H0","H0_lower","H0_upper")]
  # Describing CIs
  Hstat<-c(unlist(resH[1,]))[c(1,2,3)]
  Hstat<-round(Hstat,2)
  a<-paste0(Hstat[1]," (")
  a<-paste0(a,Hstat[2])
  a<-paste0(a,",")
  a<-paste0(a,Hstat[3])
  Hc_res1<-paste0(a,")")
  resH <- Hc_estimates[,c("H2","H2_lower","H2_upper")]
  # Describing CIs
  Hstat<-c(unlist(resH[1,]))[c(1,2,3)]
  Hstat<-round(Hstat,2)
  a<-paste0(Hstat[1]," (")
  a<-paste0(a,Hstat[2])
  a<-paste0(a,",")
  a<-paste0(a,Hstat[3])
  Hc_res2<-paste0(a,")")
  if(details){
    cat("H un-adjusted estimates-----:   ",c(H_res1), "\n")
    cat("H bias-corrected estimates--:   ",c(H_res2), "\n")
    cat("H^c un-adjusted estimates---:   ",c(Hc_res1), "\n")
    cat("H^c bias-corrected estimates:   ",c(Hc_res2), "\n")
  }

SG_CIs <- list(H_raw=H_res1,H_bc=H_res2,Hc_raw=Hc_res1,Hc_bc=Hc_res2)
 
  if(est.scale=="1/hr"){
    # Adjusted CIs
    # Recommended hr_1a (Adjusted CI for recommended)
    # Re-defined treatment.recommend flag in fs.est$df.est
    # That is, going back to estimation dataset to redefine recommendation interpretation
    fs.est$df.est[,"treat.recommend2"] <- 1 -fs.est$df.est[,"treat.recommend"]
    # treat.recommend2: 1 -- > recommend
    
    # For bootstraps above, the recommendation flag is NOT modified
    # Re-define bias-corrected estimates to correspond to treat.recommend2
    hr_1a <- SG_CIs$H_bc # recommend
    hr_0a <- SG_CIs$Hc_bc
    FSsg_tab<-SG_table(df=fs.est$df.est,SG_flag="treat.recommend2",draws=0, details=FALSE,
                    outcome.name=fs.est$outcome.name,event.name=fs.est$event.name,treat.name=fs.est$treat.name,
                    hr_1a=hr_1a, hr_0a=hr_0a, est.scale="1/hr",sg0_name="Questionable",sg1_name="Recommend")
  }
  if(est.scale=="hr") {
    hr_1a <- SG_CIs$Hc_bc
    hr_0a <- SG_CIs$H_bc
    FSsg_tab<-SG_table(df=fs.est$df.est,SG_flag="treat.recommend",draws=0, details=FALSE,
                    outcome.name=fs.est$outcome.name,event.name=fs.est$event.name,treat.name=fs.est$treat.name,
                    hr_1a=hr_1a, hr_0a=hr_0a, est.scale="hr",sg0_name="Questionable",sg1_name="Recommend")
  }  
  out <- list(results=results,SG_CIs=SG_CIs,FSsg_tab=FSsg_tab$res_out,Ystar_mat=Ystar_mat,H_estimates=H_estimates,Hc_estimates=Hc_estimates)
  }  
return(out)
}






