## Forest search cross-validation 

## Necessary packages
# Add plyr needed for rbind.fill
requiredPackages <- c("doRNG","doFuture","plyr","stringr")
## Note: cli only needed for code checking and printing when details=TRUE

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)>0){
    install.packages(new.pkg, dependencies = TRUE, quietly=TRUE, verbose=FALSE, logical.return=FALSE)
    sapply(pkg, require, character.only = TRUE, quietly=TRUE)
  }
}
ipak(requiredPackages)
sapply(c(requiredPackages), require, character.only = TRUE, quietly=TRUE)

#rbind_with_fill <- function(...) rbind(...,fill = TRUE)

forestsearch_Kfold <- function(fs.est,Kfolds=nrow(fs.est$df.est),seedit=8316951,sg0.name="Not recommend",sg1.name="Recommend",details=FALSE){
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
t.startk<-proc.time()[3]
est.scale <- c(fs.est$est.scale)
confounders.name <- c(fs.est$Allconfounders.name)
outcome.name<-c(fs.est$outcome.name)
event.name<-c(fs.est$event.name)
id.name<-c(fs.est$id.name)
treat.name<-c(fs.est$treat.name)
## Note: these names are global variables
## Extract necessary data
dfa <- fs.est$df.est[,c(confounders.name,outcome.name,event.name,id.name,treat.name,"treat.recommend")]
## re-name treat.recommend to original to compare with k-folds
names(dfa)[names(dfa) == "treat.recommend"] <- "treat.recommend.original"

dfnew=as.data.frame(dfa)

if(Kfolds==nrow(fs.est$df.est)){
df_scrambled <- copy(dfnew)
}
if(Kfolds < nrow(fs.est$df.est)){
df_scrambled <- copy(dfnew)  
set.seed(seedit)
id_sample <- sample(c(1:nrow(df_scrambled)),replace=FALSE)  
df_scrambled <- df_scrambled[id_sample,]
}

folds <- cut(seq(1,nrow(df_scrambled)),breaks=Kfolds,labels=FALSE)

if(details) cat("Range of unique left-out fold sample sizes",c(range(unique(table(folds)))),"\n")

resCV <- foreach(
  cv_index = seq_len(Kfolds),
  .options.future=list(seed=TRUE),
  .combine = "rbind.fill", 
  .errorhandling="pass"
) %dofuture% {
  testIndexes <- which(folds==cv_index,arr.ind=TRUE)
  x.test <- df_scrambled[testIndexes, ]
  x.train <- df_scrambled[-testIndexes, ]
  
 fs.train <- try(forestsearch(df.analysis=x.train, 
 is.RCT=fs.est$is.RCT, sg_focus=fs.est$sg_focus, est.scale=est.scale,
 defaultcut_names=fs.est$defaultcut_names, cut_type=fs.est$cut_type,
  Allconfounders.name=fs.est$Allconfounders.name,details=FALSE,plot.sg=FALSE,
 use_lasso=fs.est$use_lasso, use_grf=fs.est$use_grf, use_grf_only=fs.est$use_grf_only,
 dmin.grf=fs.est$dmin.grf, frac.tau=fs.est$frac.tau, conf_force=fs.est$conf_force, 
 outcome.name=fs.est$outcome.name,treat.name=fs.est$treat.name,
 event.name=fs.est$event.name,id.name=fs.est$id.name,n.min=fs.est$n.min,
 hr.threshold=fs.est$hr.threshold,hr.consistency=fs.est$hr.consistency,
 fs.splits=fs.est$fs.splits,d0.min=fs.est$d0.min,d1.min=fs.est$d1.min,
 pstop_futile=fs.est$pstop_futile,pconsistency.threshold=fs.est$pconsistency.threshold, 
 stop.threshold=fs.est$stop.threshold,max.minutes=fs.est$max.minutes,max_n_confounders=fs.est$max_n_confounders,
 replace_med_grf=fs.est$replace_med_grf,maxk=fs.est$maxk,by.risk=12),TRUE)

 if(!inherits(fs.train,"try-error") & !is.null(fs.train$sg.harm)){
    sg1 <- fs.train$sg.harm[1]
    sg2 <- fs.train$sg.harm[2]
    df.test <- get_dfpred(df.predict=x.test,sg.harm=fs.train$sg.harm,version=2)
    df.test$cvindex <- rep(cv_index,nrow(df.test))
    df.test$sg1 <- rep(sg1,nrow(df.test))
    df.test$sg2 <- rep(sg2,nrow(df.test))
    dfres <- data.table(df.test[,c(id.name,outcome.name,event.name,treat.name,"treat.recommend","treat.recommend.original","cvindex","sg1","sg2")])
  } else {
    df.test <- x.test
    df.test$cvindex <- rep(cv_index,nrow(df.test))
    df.test$sg1 <- rep(NA,nrow(df.test))
    df.test$sg2 <- rep(NA,nrow(df.test))
    df.test$treat.recommend <- rep(1.0,nrow(df.test))
    dfres <- data.table(df.test[,c(id.name,outcome.name,event.name,treat.name,"treat.recommend","treat.recommend.original","cvindex","sg1","sg2")])
 }
   return(dfres)
}  
t.now<-proc.time()[3]
t.min<-(t.now-t.startk)/60

if(length(unique(resCV$id)) != nrow(df_scrambled)) stop("K-fold cross-validation sample size does not equal full analysis dataset")
chk1 <- c(range(unique(table(folds))))
chk2 <- c(range(unique(table(resCV$cvindex))))
if(!all(chk1  == chk2)) stop("Mismatch on range of fold sample sizes")
if(details) cat("Minutes for Cross-validation: (Kfolds,minutes)=",c(Kfolds,t.min),"\n")
resCV <- as.data.frame(resCV)
if(est.scale=="1/hr"){
resCV$treat2 <- 1-resCV[,treat.name]
treat.name <- c("treat2")
sg0.name <- "Recommend"
sg1.name <- "Not recommend"
}
sg_found_count <- sum(ifelse(!is.na(resCV$sg1) | !is.na(resCV$sg2),1,0))
# Proportion of SG's found relative to number of folds (Kfolds)
propn_SG <- 100*round(c(sg_found_count/Kfolds),3)
if(details & Kfolds==nrow(fs.est$df.est)) cat("% of Kfolds where subgroup was found (valid for n-fold)",c(paste0(propn_SG,"%")),"\n")

return(list(resCV=resCV,confounders.name=confounders.name,outcome.name=outcome.name,event.name=event.name,timing_minutes=t.min,
            prop_SG_found=propn_SG,
            treat.name=treat.name,sg_analysis=fs.est$sg.harm,sg0.name=sg0.name,sg1.name=sg1.name,Kfolds=Kfolds))
}


cv_forparallel<-function(cv_index){
  stop("Need revision try error handling for fs.train")
  set.seed(8316951)
  testIndexes <- which(folds==cv_index,arr.ind=TRUE)
  x.test <- df_scrambled[testIndexes, ]
  x.train <- df_scrambled[-testIndexes, ]
 
  fs.train <- forestsearch(df.analysis=x.train, 
                           is.RCT=fs.est$is.RCT,
                           defaultcut_names=fs.est$defaultcut_names, cut_type=fs.est$cut_type,
                           sg_focus=fs.est$sg_focus,
                           Allconfounders.name=fs.est$Allconfounders.name,
                           details=FALSE,use_lasso=fs.est$use_lasso, use_grf=fs.est$use_grf, use_grf_only=fs.est$use_grf_only,
                           dmin.grf=fs.est$dmin.grf, frac.tau=fs.est$frac.tau,
                           conf_force=fs.est$conf_force, outcome.name=fs.est$outcome.name,treat.name=fs.est$treat.name,
                           event.name=fs.est$event.name,id.name=fs.est$id.name,
                           n.min=fs.est$n.min,hr.threshold=fs.est$hr.threshold,hr.consistency=fs.est$hr.consistency,
                           fs.splits=fs.est$fs.splits,d0.min=fs.est$d0.min,d1.min=fs.est$d1.min,
                           pstop_futile=fs.est$pstop_futile,pconsistency.threshold=fs.est$pconsistency.threshold, 
                           stop.threshold=fs.est$stop.threshold,
                           max.minutes=fs.est$max.minutes,
                           est.scale=est.scale,
                           max_n_confounders=fs.est$max_n_confounders,
                           replace_med_grf=fs.est$replace_med_grf,
                           maxk=fs.est$maxk,by.risk=12,
                           plot.sg=FALSE)
  
  if(!is.null(fs.train$sg.harm)){
    sg1 <- fs.train$sg.harm[1]
    sg2 <- fs.train$sg.harm[2]
    df.test <- get_dfpred(df.predict=x.test,sg.harm=fs.train$sg.harm,version=2)
    df.test$cvindex <- rep(cv_index,nrow(df.test))
    df.test$sg1 <- rep(sg1,nrow(df.test))
    df.test$sg2 <- rep(sg2,nrow(df.test))
    dfres <- data.table(df.test[,c(id.name,outcome.name,event.name,treat.name,"treat.recommend","treat.recommend.original","cvindex","sg1","sg2")])
  }
  
  if(is.null(fs.train$sg.harm)){
    df.test <- x.test
    df.test$cvindex <- rep(cv_index,nrow(df.test))
    df.test$sg1 <- rep(NA,nrow(df.test))
    df.test$sg2 <- rep(NA,nrow(df.test))
    df.test$treat.recommend <- rep(1.0,nrow(df.test))
    dfres <- data.table(df.test[,c(id.name,outcome.name,event.name,treat.name,"treat.recommend","treat.recommend.original","cvindex","sg1","sg2")])
  }
  return(dfres)
}  



forestsearch_KfoldOut<-function(res,details=FALSE,outall=FALSE){
  
confounders.name <- c(res$confounders.name)
outcome.name <- c(res$outcome.name)
event.name <- c(res$event.name)
treat.name <- c(res$treat.name)
sg_analysis <- c(res$sg_analysis)
sg0.name <- c(res$sg0.name)
sg1.name <- c(res$sg1.name)

Kfolds <- res$Kfolds
df_CV <- res$resCV

## outall=FALSE --> only output find_metrics and sens_metrics_original
## Extract sg1 and sg2
  sg1 <- sg2 <- rep(NA,Kfolds)
  for(ks in 1:Kfolds){
    ## First element since all identical for same cvindex (=ks)  
    sg1[ks] <- subset(df_CV, cvindex==ks)$sg1[1] 
    sg2[ks] <- subset(df_CV, cvindex==ks)$sg2[1] 
  }
  SGs_found <- cbind(sg1,sg2)

  CV_summary <- CV_sgs(sg1=sg1, sg2=sg2, confs=confounders.name, sg_analysis=sg_analysis)
  
  if(details){
    cat("Any found",c(mean(CV_summary$any_found)),"\n")
    cat("Exact match",c(mean(CV_summary$exact_match)),"\n")
    cat("At least 1 match",c(mean(CV_summary$one_match)),"\n")
    cat("Cov 1 any",c(mean(CV_summary$cov1_any)),"\n")
    cat("Cov 2 any",c(mean(CV_summary$cov2_any)),"\n")
    cat("Cov 1 and 2 any",c(mean(ifelse(CV_summary$cov1_any & CV_summary$cov2_any,1,0))),"\n")
    cat("Cov 1 exact",c(mean(CV_summary$cov1_exact)),"\n")
    cat("Cov 2 exact",c(mean(CV_summary$cov2_exact)),"\n")
  }
  
  find_metrics <- c(mean(CV_summary$any_found),mean(CV_summary$exact_match),mean(CV_summary$one_match),
                    mean(CV_summary$cov1_any),mean(CV_summary$cov2_any),mean(ifelse(CV_summary$cov1_any & CV_summary$cov2_any,1,0)),
                    mean(CV_summary$cov1_exact),mean(CV_summary$cov2_exact))
  
  names(find_metrics) <- c("Any","Exact","At least 1","Cov1","Cov2","Cov 1 & 2","Cov1 exact","Cov2 exact")
  
  # Propn agreement in H and H^c
  # between sample estimate (original) and cross-validation
  # Possible all folds return no subgroup
  n_sgfound <- length(unique(df_CV$treat.recommend))
  if(n_sgfound==2){
  tabit <- with(df_CV,table(treat.recommend,treat.recommend.original))
  sensH <- tabit[1,1]/sum(tabit[,1])
  sensHc <- tabit[2,2]/sum(tabit[,2])
  ppvH <- tabit[1,1]/sum(tabit[1,])
  ppvHc <- tabit[2,2]/sum(tabit[2,])
  }
  if(n_sgfound==1){
    tabit <- with(df_CV,table(treat.recommend,treat.recommend.original))
    sensH <- tabit[1,1]/sum(tabit[,1])
    sensHc <- NA
    ppvH <- tabit[1,1]/sum(tabit[1,])
    ppvHc <- NA
  }
  sens_metrics_original <- c(sensH,sensHc,ppvH,ppvHc)
  names(sens_metrics_original) <- c("sens_H","sens_Hc","ppv_H","ppv_Hc")
  
  if(details) cat("Agreement (sens, ppv) in H and Hc:",c(sensH,sensHc,ppvH,ppvHc),"\n")
  
  if(outall){
  itt_tab <- SG_table(df=as.data.frame(df_CV),SG_flag="ITT",outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,draws=0)
  if(n_sgfound==2){
  SG_tab_Kfold <- SG_table(df=as.data.frame(df_CV),SG_flag="treat.recommend",sg1_name=sg1.name,sg0_name=sg0.name,
  outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,draws=0)
  }
  if(n_sgfound==1){
  SG_tab_Kfold <- itt_tab  
  }
 # Data analysis subgroups
 SG_tab_original <- SG_table(df=as.data.frame(df_CV),SG_flag="treat.recommend.original",sg1_name=sg1.name,sg0_name=sg0.name,
 outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,draws=0)

 # Combine tables starting with ITT
if(n_sgfound==2){ 
temp1 <- itt_tab$res_out[,c("Subgroup","n","n1","n0","m1","m0","RMST","HR")]
temp2a <- SG_tab_original$res_out[1,c("Subgroup","n","n1","n0","m1","m0","RMST","HR")]
temp2b <- SG_tab_original$res_out[2,c("Subgroup","n","n1","n0","m1","m0","RMST","HR")]
temp3a <- SG_tab_Kfold$res_out[1,c("Subgroup","n","n1","n0","m1","m0","RMST","HR")]
temp3b <- SG_tab_Kfold$res_out[2,c("Subgroup","n","n1","n0","m1","m0","RMST","HR")]
tab_all <- rbind(temp1,temp2a,temp3a,temp2b,temp3b)
colnames(tab_all) <- c("Subgroup","n","n1","n0","m1","m0","RMST","Hazard ratio") 
rownames(tab_all) <- c("Overall","FA_0","KfA_0","FA_1","KfA_1")
}
 
 if(n_sgfound==1){ 
   temp1 <- itt_tab$res_out[,c("Subgroup","n","n1","n0","m1","m0","RMST","HR")]
   temp2a <- SG_tab_original$res_out[1,c("Subgroup","n","n1","n0","m1","m0","RMST","HR")]
   temp2b <- SG_tab_original$res_out[2,c("Subgroup","n","n1","n0","m1","m0","RMST","HR")]
   tab_all <- rbind(temp1,temp2a,temp2b)
   colnames(tab_all) <- c("Subgroup","n","n1","n0","m1","m0","RMST","Hazard ratio") 
   rownames(tab_all) <- c("Overall","FA_0","FA_1")
 }
 

if(details) print(tab_all)

out <- list(itt_tab=itt_tab,SG_tab_original=SG_tab_original,SG_tab_Kfold=SG_tab_Kfold,
CV_summary=CV_summary,sens_metrics_original=sens_metrics_original,find_metrics=find_metrics,
SGs_found=SGs_found,tab_all=tab_all)
}
 if(!outall){ 
   out <- list(sens_metrics_original=sens_metrics_original,find_metrics=find_metrics)
 }
 return(out)
  }
  
## Summary of matches with full data analysis subgroups

CV_sgs <- function(sg1,sg2,confs,sg_analysis){
  
any_found <- ifelse(!is.na(sg1) | !is.na(sg2),1,0)
sg_depth <- length(sg_analysis)

if(sg_depth==2){
  sg1a <- sg_analysis[1]
  sg2a <- sg_analysis[2]
  ## Exact match on both to analysis data
  exact_match <- ifelse((sg1 == sg1a | sg2 == sg1a) & (sg1 == sg2a | sg2 == sg2a),1,0)
  exact_match[is.na(exact_match)] <- 0.0
  ## At least 1 exact match
  one_match <- ifelse((sg1 == sg1a | sg2 == sg1a) | (sg1 == sg2a | sg2 == sg2a),1,0)
  one_match[is.na(one_match)] <- 0.0
  ## Cov 1 exact
  cov1_match <- ifelse((sg1 == sg1a | sg2 == sg1a),1,0)
  cov1_match[is.na(cov1_match)] <- 0.0
  cov2_match <- ifelse((sg1 == sg2a | sg2 == sg2a),1,0)
  cov2_match[is.na(cov2_match)] <- 0.0
  ## Find confounder names involved in sg1a and sg2a
  ## Add { or !{ to names for matching (a bit tedious, but let's see)
  dda <- charmatch("{",sg1a, nomatch=0)
  ddb <- charmatch("!{",sg1a, nomatch=0)
  if(dda ==1) aa <- rep("{",length(confs))
  if(ddb ==1) aa <- rep("!{",length(confs))
  
  temp <- paste0(aa,confs)
  temp <-paste0(temp,"}")
  ## Is sg1a confounder (NOT necessarily same cut) involved in any
  loc_name <- charmatch(temp,sg1a)
  index_name <- which(loc_name==1)

  if(length(index_name)==0){
  temp <- paste0(aa,confs)
  loc_name <- charmatch(temp,sg1a)
  index_name <- which(loc_name==1)
  }
  
  # Find exact match
  if(length(index_name)>1){
    confs2 <- confs[index_name]
    lc <- str_length(confs2)
   if(dda==1) ctoget <- str_sub(sg1a,2,max(lc))
   if(ddb==1) ctoget <- str_sub(sg1a,3,max(lc))
    itoget <- which(confs == ctoget)
    cfs <- confs[itoget]
  }
  if(length(index_name)==1) cfs <- confs[which(loc_name==1)]

  bb1 <- grepl(cfs,sg1)
  bb2 <- grepl(cfs,sg2)
  cov1_any <- ifelse(bb1 | bb2, 1,0)
  rm("bb1","bb2","cfs")

   
  ## Second
  dda <- charmatch("{",sg2a, nomatch=0)
  ddb <- charmatch("!{",sg2a, nomatch=0)
  if(dda ==1) aa <- rep("{",length(confs))
  if(ddb ==1) aa <- rep("!{",length(confs))
  temp <- paste0(aa,confs)
  temp <-paste0(temp,"}")
  
  loc_name <- charmatch(temp,sg2a)
  index_name <- which(loc_name==1)
  
  if(length(index_name)==0){
    temp <- paste0(aa,confs)
    loc_name <- charmatch(temp,sg2a)
    index_name <- which(loc_name==1)
  }
  
  # Find exact match
  if(length(index_name)>1){
    confs2 <- confs[index_name]
    lc <- str_length(confs2)
    if(dda==1) ctoget <- str_sub(sg2a,2,max(lc))
    if(ddb==1) ctoget <- str_sub(sg2a,3,max(lc))
    itoget <- which(confs == ctoget)
    cfs <- confs[itoget]
  }
  if(length(index_name)==1) cfs <- confs[which(loc_name==1)]

  bb1 <- grepl(cfs,sg1)
  bb2 <- grepl(cfs,sg2)
  cov2_any <- ifelse(bb1 | bb2, 1,0)
}
if(sg_depth==1){
  sg1a <- sg_analysis[1]
  ## Exact match on both to analysis data
  exact_match <- ifelse((sg1 == sg1a | sg2 == sg1a),1,0)
  exact_match[is.na(exact_match)] <- 0.0
  dda <- charmatch("{",sg1a, nomatch=0)
  ddb <- charmatch("!{",sg1a, nomatch=0)
  if(dda ==1) aa <- rep("{",length(confs))
  if(ddb ==1) aa <- rep("!{",length(confs))
  temp <- paste0(aa,confs)
  temp <-paste0(temp,"}")
  ## Is sg1a confounder (NOT necessarily same cut) involved in any
  loc_name <- charmatch(temp,sg1a)
  index_name <- which(loc_name==1)

  if(length(index_name)==0){
    temp <- paste0(aa,confs)
    loc_name <- charmatch(temp,sg1a)
    index_name <- which(loc_name==1)
  }
  
  # If names have numbers which may not be unique
  # E.g., "z1", "z11", this will match both
  # Find exact
  # Find exact match
  if(length(index_name)>1){
  confs2 <- confs[index_name]
  lc <- str_length(confs2)
  if(dda==1) ctoget <- str_sub(sg1a,2,max(lc))
  if(ddb==1) ctoget <- str_sub(sg1a,3,max(lc))
  itoget <- which(confs == ctoget)
  cfs <- confs[itoget]
  }
  if(length(index_name)==1) cfs <- confs[which(loc_name==1)]
  bb1 <- grepl(cfs,sg1)
  bb2 <- grepl(cfs,sg2)
  cov1_any <- ifelse(bb1 | bb2, 1,0)
  one_match <- exact_match
  cov2_any <- NA

  cov1_match <- exact_match
  cov2_match <- NA
  
}
return(list(any_found=any_found, exact_match=exact_match, one_match=one_match, cov1_any=cov1_any, cov2_any=cov2_any, cov1_exact=cov1_match, cov2_exact=cov2_match))
}

# Simulation 10-fold cross-validation

forestsearch_tenfold <- function(fs.est,sims,Kfolds=10,details=TRUE){
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
t.start<-proc.time()[3]
## Will summarize cross-validated metrics
## Initialize here (will be concatenated below)
sens_out <- NULL
find_out <- NULL
confounders.name <- c(fs.est$Allconfounders.name)
outcome.name<-c(fs.est$outcome.name)
event.name<-c(fs.est$event.name)
id.name<-c(fs.est$id.name)
treat.name<-c(fs.est$treat.name)
sg_analysis <- c(fs.est$sg.harm)
dfa <- fs.est$df.est[,c(confounders.name,outcome.name,event.name,id.name,treat.name,"treat.recommend")]
## re-name treat.recommend to original to compare with k-folds
names(dfa)[names(dfa) == "treat.recommend"] <- "treat.recommend.original"
dfnew <- as.data.frame(dfa)
for(ksim in 1:sims){
  t.startk<-proc.time()[3]
  df_scrambled = copy(dfnew)
  set.seed(8316951+1000*ksim)
  id_sample <- sample(c(1:nrow(df_scrambled)),replace=FALSE)
  df_scrambled=df_scrambled[id_sample,]
  folds=cut(seq(1,nrow(df_scrambled)),breaks=Kfolds,labels=FALSE)
  if(details & ksim<=3) cat("Range of unique left-out fold sample sizes (first 3 simulations)",c(range(unique(table(folds)))),"\n")
   resCV <- foreach(
    cv_index = seq_len(Kfolds),
    .options.future=list(seed=TRUE,add=c("df_scrambled","folds","fs.est")),
    .combine=rbind.fill,
    .errorhandling="pass"
  ) %dofuture% {
    testIndexes <- which(folds==cv_index,arr.ind=TRUE)
    x.test <- df_scrambled[testIndexes, ]
    x.train <- df_scrambled[-testIndexes, ]
    
    fs.train <- try(forestsearch(df.analysis=x.train, details=FALSE,
                             is.RCT=fs.est$is.RCT, 
                             defaultcut_names=fs.est$defaultcut_names, cut_type=fs.est$cut_type,
                             sg_focus=fs.est$sg_focus,
                             Allconfounders.name=fs.est$Allconfounders.name,
                             use_lasso=fs.est$use_lasso, use_grf=fs.est$use_grf, use_grf_only=fs.est$use_grf_only,
                             dmin.grf=fs.est$dmin.grf, frac.tau=fs.est$frac.tau,
                             conf_force=fs.est$conf_force, outcome.name=fs.est$outcome.name,treat.name=fs.est$treat.name,
                             event.name=fs.est$event.name,id.name=fs.est$id.name,
                             n.min=fs.est$n.min,hr.threshold=fs.est$hr.threshold,hr.consistency=fs.est$hr.consistency,
                             fs.splits=fs.est$fs.splits,d0.min=fs.est$d0.min,d1.min=fs.est$d1.min,
                             pstop_futile=fs.est$pstop_futile,pconsistency.threshold=fs.est$pconsistency.threshold, 
                             stop.threshold=fs.est$stop.threshold,
                             max.minutes=fs.est$max.minutes,
                             max_n_confounders=fs.est$max_n_confounders,
                             replace_med_grf=fs.est$replace_med_grf,
                             est.scale=fs.est$est.scale,
                             maxk=fs.est$maxk,by.risk=12,
                             plot.sg=FALSE),TRUE)
    
if(!inherits(fs.train,"try-error")){    
    if(!is.null(fs.train$sg.harm)){
     if(details & ksim <= 3){
     cat("Simulation, Fold =",c(ksim,cv_index),"\n")  
     print(fs.train$sg.harm)  
      }      
      sg1 <- fs.train$sg.harm[1]
      sg2 <- fs.train$sg.harm[2]
      
      df.test <- get_dfpred(df.predict=x.test,sg.harm=fs.train$sg.harm,version=2)
      df.test$cvindex <- rep(cv_index,nrow(df.test))
      df.test$sg1 <- rep(sg1,nrow(df.test))
      df.test$sg2 <- rep(sg2,nrow(df.test))
      dfres <- data.table(df.test[,c(id.name,outcome.name,event.name,treat.name,"treat.recommend","treat.recommend.original","cvindex","sg1","sg2")])
     }  else {
     if(details & ksim <= 3) cat("Subgroup not found (return ITT), Fold =",c(ksim,cv_index),"\n")  
      df.test <- x.test
      df.test$cvindex <- rep(cv_index,nrow(df.test))
      df.test$sg1 <- rep(NA,nrow(df.test))
      df.test$sg2 <- rep(NA,nrow(df.test))
      df.test$treat.recommend <- rep(1.0,nrow(df.test))
      dfres <- data.table(df.test[,c(id.name,outcome.name,event.name,treat.name,"treat.recommend","treat.recommend.original","cvindex","sg1","sg2")])
    }
}
if(inherits(fs.train,"try-error")){    
  # cat("Hey look at me","\n")
  # cat("Simulation, Fold =",c(ksim,cv_index),"\n")  
  
  df.test <- x.test
  df.test$cvindex <- rep(cv_index,nrow(df.test))
  df.test$sg1 <- rep(NA,nrow(df.test))
  df.test$sg2 <- rep(NA,nrow(df.test))
  df.test$treat.recommend <- rep(1.0,nrow(df.test))
  dfres <- data.table(df.test[,c(id.name,outcome.name,event.name,treat.name,"treat.recommend","treat.recommend.original","cvindex","sg1","sg2")])
}  
  return(dfres)
  }  
  if(length(unique(resCV$id)) != nrow(df_scrambled)) stop("K-fold cross-validation sample size does not equal full analysis dataset")
  chk1 <- c(range(unique(table(folds))))
  chk2 <- c(range(unique(table(resCV$cvindex))))
  if(!all(chk1  == chk2)) stop("Mismatch on range of fold sample sizes")
  
  res <- list(resCV=resCV,Kfolds=Kfolds,confounders.name=confounders.name,
  outcome.name=outcome.name,event.name=event.name,id.name=id.name,
  treat.name=treat.name,sg_analysis=sg_analysis,sg0.name="Not recommend",sg1.name="Recommend")
  # Note: sg0.name and sg1.name not used here
  out <- forestsearch_KfoldOut(res=res,outall=FALSE,details=FALSE)
  ## Projected timing
   if(details & (ksim %in% c(c(1:5),10,20,50,100,200,300))){
     cat("Kfold iteration done=",c(ksim),"\n")
   }
 sens_out <- rbind(sens_out,out$sens_metrics_original)
 find_out <- rbind(find_out,out$find_metrics)
}
t.now<-proc.time()[3]
t.min<-(t.now-t.start)/60
if(details){
cat("Minutes for Cross-validation: Kfolds, sims, minutes=",c(Kfolds,sims,t.min),"\n")
cat("Projected hours per 100 sims",c((t.min/60)*(100/sims)),"\n")
}
sens_summary <- apply(sens_out,2,median,na.rm=TRUE)
find_summary <- apply(find_out,2,median,na.rm=TRUE)
return(list(sens_summary=sens_summary,find_summary=find_summary,sens_out=sens_out,find_out=find_out))
}


