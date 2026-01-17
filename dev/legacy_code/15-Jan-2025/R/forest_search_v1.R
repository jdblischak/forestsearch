
# Necessary packages
# 8-28-2024 remove plyr and dplyr

# Removed restrictions: (19-Sept-2024)
#if(any(grep("_",confounders.name)) | any(grep("-",confounders.name))) stop("confounders names with underscores or hyphenation not allowed (please change [eg] var_1 to var1")  
#if(any(grep("_",conf_force)) | any(grep("-",conf_force))) stop("confounders/conf_force names with underscores or hyphenation not allowed (please change [eg] var_1 to var1")  

# Replace confounders.name with confounders.name (17-March-2025)



requiredPackages <- c("grf","policytree","data.table","randomForest","survival","glmnet")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)>0){
    install.packages(new.pkg, dependencies = TRUE, quietly=TRUE, verbose=FALSE, logical.return=FALSE)
  sapply(pkg, require, character.only = TRUE, quietly=TRUE)
  }
}
ipak(requiredPackages)
sapply(c(requiredPackages), require, character.only = TRUE, quietly=TRUE)

# Prediction dataset
get_dfpred<-function(df.predict,sg.harm,version=1){
  # Prediction dataset for cross-validation:
  # df.train is training dataset for estimation, df.pred=df.test 
  # Or just original analysis dataset which facilitates summaries based on 
  # source factors and comparison with other methods  
if(version==1)  df.pred<-dummy(df.predict)
if(version==2) df.pred<-dummy2(df.predict)

  df.pred$treat.recommend<-NA
  id.harm<-c(paste(sg.harm, collapse="==1 & "))
  id.harm<-paste(id.harm,"==1")
  df.pred.0<-subset(df.pred,eval(parse(text=id.harm)))
  
  if(nrow(df.pred.0)>0){
  df.pred.0$treat.recommend<-0
  }
  id.noharm<-c(paste(sg.harm, collapse="!=1 | "))
  id.noharm<-paste(id.noharm,"!=1")
  df.pred.1<-subset(df.pred,eval(parse(text=id.noharm)))
  if(nrow(df.pred.1)>0){
    df.pred.1$treat.recommend<-1
  }
if(nrow(df.pred.0)>0 & nrow(df.pred.1)>0)  df.pred.out<-data.frame(rbind(df.pred.0,df.pred.1))
if(nrow(df.pred.0)==0 & nrow(df.pred.1)>0)  df.pred.out<-data.frame(df.pred.1)
if(nrow(df.pred.0)>0 & nrow(df.pred.1)==0)  df.pred.out<-data.frame(df.pred.0)
return(df.pred.out)
}


forestsearch<-function(df.analysis=NULL,df.predict=NULL,df.test=NULL,
is.RCT=TRUE,seedit=8316951,
est.scale="hr",
confounders.name=NULL,
use_lasso=FALSE,
use_grf=TRUE, use_grf_only=FALSE,
grf_res=NULL,grf_cuts=NULL,max_n_confounders=1000,
grf_depth=2, dmin.grf=12, frac.tau=0.6,
conf_force=NULL,
defaultcut_names=NULL, cut_type="default",
replace_med_grf=FALSE,
cont.cutoff=4,
conf.cont_medians=NULL,
conf.cont_medians_force=NULL,
outcome.name=NULL,event.name=NULL,treat.name=NULL,id.name=NULL,
n.min=60,hr.threshold=1.25,hr.consistency=1.0,
sg_focus="hr",
fs.splits=1000,m1.threshold=Inf,
stop.threshold=0.90,pconsistency.threshold=0.90,
showten_subgroups=FALSE,
d0.min=10,d1.min=10,
pstop_futile=0.5,max.minutes=3,minp=0.025,
details=FALSE,maxk=2,by.risk=3,plot.sg=TRUE,
vi.grf.min=-0.2){

if(showten_subgroups){
# reset stop.threshold (>1) to override stopping
stop.threshold <- 1.1
}  
  
if(!(sg_focus %in% c("hr","Nsg", "Msg", "Nsg_only", "Msg_only"))) stop("sg_focus must be either hr, Nsg (Nsg_only), or Msg (Msg_only)")

if(!(cut_type %in% c("default","median"))) stop("only default and median cuts")  
  
if(all(confounders.name %in% names(df.analysis)) != TRUE) stop("Not all confounders found in dataset")   

if(!is.null(defaultcut_names)){  
if(all(defaultcut_names %in% names(df.analysis)) != TRUE) stop("Not all confounders for default cuts found in dataset")   
}

if(is.null(id.name)){  
# Add id 
df.analysis$id <- c(1:nrow(df.analysis))
id.name <- c("id")
}  

if(is.null(df.analysis)) stop("Dataset (df.analysis) required")  
if(is.null(confounders.name)) stop("Confounder names (confounders.name) required")
if(is.null(outcome.name) | is.null(event.name) | is.null(treat.name)) stop("outcome.name, event.name, and treat.name required (missing at least 1)")
if(is.null(hr.threshold) | is.null(hr.consistency) | is.null(pconsistency.threshold)) stop("hr.threshold, hr.consistency, pconsistency.threshold required (missing at least 1)")
if(is.null(stop.threshold) | is.null(pstop_futile)) stop("stop.threshold, pstop_futile required)")  
  
if(sg_focus=="Nsg_only" & stop.threshold > pconsistency.threshold){
stop.threshold <- pconsistency.threshold
warning("For Nsg_only stop.threshold is reset to equal pconsistency.threshold")
}  
if(sg_focus=="Nsg" & stop.threshold < 0.99){
stop.threshold <- 0.99
warning("For NSG re-setting stop.threshold to 0.99")
}  
  
  
# Create id by observation number if id.name is not provided  
if(is.null(id.name)){
df.analysis$id <- c(1:nrow(df.analysis))
id.name <- c("id")
}
# Sort data by id
df.analysis <- df.analysis[order(df.analysis[,id.name]),]
# Check and omit any misssing data
var_names<-c(confounders.name,outcome.name,event.name,id.name,treat.name)
dfa <- df.analysis %>% select(one_of(var_names))
# How many subjects excluded if omitting missing data
df.analysis <- na.omit(dfa)
if(nrow(df.analysis) != nrow(dfa)){
cat("Total excluded by omitting missing data=",c(nrow(dfa)-nrow(na.omit(dfa))),"\n")
}

t.start_all<-proc.time()[3]
# Set use_grf if use_grf_only and set lasso to FALSE
if(use_grf_only){
use_grf <- TRUE
use_lasso <- FALSE
}

# If grf_res is populated
if(use_grf & !is.null(grf_res)){
if(is.null(grf_res$tree.cuts)){
# If no cuts, then re-set use_grf
use_grf <- FALSE
}  
if(!is.null(grf_res$tree.cuts)){
# GRF cuts
grf_cuts <- grf_res$tree.cuts
}  
}
  
# If using grf and not populated then run grf
if(use_grf & is.null(grf_res)){
if(details){
cat("FS: GRF stage for cut selection with dmin,tau=",c(dmin.grf,frac.tau),"\n")
}
grf_res <- NULL  
grf_res<-try(grf.subg.harm.survival(data=df.analysis,confounders.name=confounders.name,
outcome.name=outcome.name,RCT=is.RCT,seedit=seedit,maxdepth=grf_depth,
event.name=event.name,id.name=id.name,treat.name=treat.name,n.min=n.min,dmin.grf=dmin.grf,
frac.tau=frac.tau,details=details),TRUE)

# If grf returned but NO sg found set use_grf=FALSE
if(!inherits(grf_res,"try-error") & is.null(grf_res$sg.harm)){
use_grf <- FALSE
if(details) cat("NO GRF cuts meeting delta(RMST): dmin.grf=",c(dmin.grf),"\n")
}
# If grf returned and sg found 
if(!inherits(grf_res,"try-error") & !is.null(grf_res$sg.harm)){
  
plot(grf_res$tree,leaf.labels=c("Control","Treat"))

grf_cuts <- grf_res$tree.cuts
}
}

FSdata <-try(get_FSdata(df=df.analysis,confounders.name=confounders.name,
outcome.name=outcome.name, event.name=event.name,
use_lasso=use_lasso,use_grf=use_grf, use_grf_only=use_grf_only,
grf_cuts=grf_cuts,
cut_type=cut_type, defaultcut_names=defaultcut_names,
replace_med_grf=replace_med_grf,
cont.cutoff=cont.cutoff,                     
conf.cont_medians=conf.cont_medians,conf.cont_medians_force=conf.cont_medians_force,
conf_force=conf_force,
details=details),TRUE)

if(inherits(FSdata,"try-error")) stop("FSdata error")  
  
if(!inherits(FSdata,"try-error")){  
lassoomit <- FSdata$lassoomit
lassokeep <- FSdata$lassokeep
df <- FSdata$df
FSconfounders.name <- FSdata$confs_names
confs_labels <- FSdata$confs
if(is.null(df.predict)) df.predict <- df
  if(!is.null(vi.grf.min)){
  # Use GRF for screening and ordering
  # Covariates need to be converted to numeric scale
  X<-as.matrix(df[,FSconfounders.name])
  # Convert to numeric
  X<-apply(X,2,as.numeric)
  Y<-df[,outcome.name]
  W<-df[,treat.name]
  D<-df[,event.name]
  tau.rmst<-min(c(max(Y[W==1 & D==1]),max(Y[W==0 & D==1])))
  # For screening we take 0.9*tau.rms
  if(!is.RCT) cs.forest<-try(suppressWarnings(causal_survival_forest(X,Y,W,D,horizon=0.9*tau.rmst,seed=8316951)),TRUE)
  if(is.RCT) cs.forest<-try(suppressWarnings(causal_survival_forest(X,Y,W,D,W.hat=0.5,horizon=0.9*tau.rmst,seed=8316951)),TRUE)
  
  vi.cs<-round(variable_importance(cs.forest),4)
  vi.cs2<-data.frame(confs_labels,FSconfounders.name,vi.cs)
  vi.order<-order(vi.cs,decreasing=TRUE)
  vi.cs2<-vi.cs2[vi.order,]
  
  conf.screen<-vi.cs2[,2]
  vi_ratio <- vi.cs2[,3] / max(vi.cs2[,3])
  selected.vars <- which(vi_ratio > vi.grf.min)
  conf.screen<-conf.screen[selected.vars]
  # Restrict to max of max_n_confounders
  # Keeping 1st lmax ordered per vi
  lmax <- min(c(length(conf.screen),max_n_confounders))
  conf.screen_limit <- conf.screen[c(1:lmax)]
  conf.screen <- conf.screen_limit
  if(details){
  cat("Number of factors evaluated=",c(lmax),"\n")
  cat("Confounders per grf screening",conf.screen,"\n")
  vi_res <- vi.cs2[selected.vars,]
  # Re-name for printing
  names(vi_res)<-c("Factors","Labels","VI(grf)")
  print(vi_res)
  }
   # Returned as confounders.evaluated
  }
  # Otherwise, do not screen and use initial factors
  if(is.null(vi.grf.min)) conf.screen<-FSconfounders.name
  df.confounders<-df[,conf.screen]
  df.confounders<-dummy(df.confounders)
  Event<-D
  Treat<-W
  # name identification as "id" for merging (df.predict) sg membership flags
  id<-df[,c(id.name)]
  df.fs<-data.frame(Y,Event,Treat,id,df.confounders)
  Z<-as.matrix(df.confounders)
  colnames(Z)<-names(df.confounders)
  
  find.grps <- subgroup.search(Y=Y,Event=Event,Treat=Treat,ID=id,Z=Z,d0.min=d0.min,d1.min=d1.min,n.min=n.min,
  hr.threshold=hr.threshold,max.minutes=max.minutes,details=details,maxk=maxk)
  
  sg.harm<-NULL 
  df.est_out<-NULL
  df.predict_out<-NULL 
  df.test_out <- NULL
  grp.consistency<-NULL
  
  max_sg_est<-find.grps$max_sg_est
  
  prop_maxk <- find.grps$prop_max_count
  
  if((!is.null(find.grps$out.found) & any(find.grps$out.found$hr.subgroups$HR>hr.consistency)) | any(find.grps$out.found$hr.subgroups$HR>hr.consistency)){# Found something?
    if(plot.sg & is.null(by.risk)) by.risk<-round(max(Y)/12,0)
    if(details) cat("# of candidate subgroups (meeting HR criteria) = ",c(nrow(find.grps$out.found$hr.subgroups)),"\n")
    
    grp.consistency<-subgroup.consistency(df=df.fs,
    hr.subgroups=find.grps$out.found$hr.subgroups,
    Lsg=find.grps$L,
    sg_focus=sg_focus,
    confs_labels=confs_labels,
    hr.threshold=hr.threshold,hr.consistency=hr.consistency,m1.threshold=m1.threshold,
    pstop_futile=pstop_futile, pconsistency.threshold=pconsistency.threshold,
    n.splits=fs.splits,maxk=maxk,
    stop.threshold=stop.threshold,details=details,by.risk=by.risk,plot.sg=plot.sg)
    
    
    if(details){
    t.end_all<-proc.time()[3]
    t.min_all<-(t.end_all-t.start_all)/60
    cat("Minutes forestsearch overall=",c(t.min_all),"\n")
  }
# Found something and a subgroup    
if(!is.null(grp.consistency$sg.harm)){   
sg.harm <- grp.consistency$sg.harm  
# data containing id and treatment flag  
temp <- grp.consistency$df_flag  
# Merge to analysis data and add treatment flag (all.x=TRUE)
df.est_out <- merge(df, temp, by="id", all.x=TRUE)  
 # Return df.predict    
if(!is.null(df.predict)){
# This does not work if df.predict is test sample in which case
# cannot match to df_flag by id
df.predict_out <- merge(df.predict, temp, by="id", all.x=TRUE)
}
if(!is.null(df.test)){
  df.test_out <- get_dfpred(df.predict=df.test,sg.harm=grp.consistency$sg.harm,version=2)
  if(est.scale=="1/hr"){
    df.test_out$treat.recommend <- 1-df.test_out$treat.recommend
}
   }
} 
 } # Found something 
out <- list(grp.consistency=grp.consistency,find.grps=find.grps,
confounders.name=confounders.name,
outcome.name=outcome.name,event.name=event.name,
treat.name=treat.name,id.name=id.name,
is.RCT=is.RCT,
max_n_confounders=max_n_confounders,
confounders.candidate=FSconfounders.name,
confounders.evaluated=confs_labels,
df.est=df.est_out,df.predict=df.predict_out,df.test=df.test_out,
grf_res=grf_res,
sg_focus=sg_focus,
sg.harm=sg.harm,
grf_cuts=grf_cuts,
conf_force=conf_force, defaultcut_names=defaultcut_names,
cut_type=cut_type,
replace_med_grf=replace_med_grf,
grf_depth=grf_depth,
dmin.grf=dmin.grf, 
frac.tau=frac.tau,
use_lasso=use_lasso, use_grf=use_grf, use_grf_only=use_grf_only,
hr.threshold=hr.threshold,
hr.consistency=hr.consistency,m1.threshold=m1.threshold,
pstop_futile=pstop_futile, 
pconsistency.threshold=pconsistency.threshold,
fs.splits=fs.splits,
d0.min=d0.min,
d1.min=d1.min,
maxk=maxk,
stop.threshold=stop.threshold,
n.min=n.min,
vi.grf.min=vi.grf.min,
frac.tau=frac.tau,
lassokeep=lassokeep,
lassoomit=lassoomit,
max.minutes=max.minutes,
prop_maxk=prop_maxk, est.scale=est.scale,
max_sg_est=max_sg_est,pstop_futile=pstop_futile)
}
# Fsdata NOT successful
if(inherits(FSdata,"try-error")){
  out <- list(sg.harm=NULL)  
}  
return(out)
}

fs.estimates.out<-function(df,fs.res=NULL,dgm=NULL,cox.formula.sim=NULL,cox.formula.adj.sim=NULL,analysis="FS",frac.tau_fs=1){
p.cens<-mean(1-df[,event.name])  

Y<-df[,outcome.name]
W<-df[,treat.name]
D<-df[,event.name]
taumax<-frac.tau_fs*min(c(max(Y[W==1 & D==1]),max(Y[W==0 & D==1])))

if(!is.null(fs.res) & !is.null(fs.res$sg.harm_label)){
sg.harm.fs<-fs.res$sg.harm_label
temp <- fs.res$df_flag  
# Merge to analysis data and add treatment flag (all.x=TRUE)
df.pred <- merge(df, temp, by="id", all.x=TRUE)  
}

if(is.null(fs.res)){
sg.harm.fs<-NULL
}

fit<-summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
hr.itt<-c(fit[1])
l.itt<-c(fit[3])
u.itt<-c(fit[4])
rm("fit")

fit<-summary(coxph(cox.formula.adj.sim,data=df,robust=FALSE))$conf.int
hr.adj.itt<-c(fit[1,1])
l.adj.itt<-c(fit[1,3])
u.adj.itt<-c(fit[1,4])
rm("fit")

if(dgm$model=="null" & !is.null(dgm$fs.harm.true)) stop("For dgm model null fs.harm.true should be null")

# dgm$fs.harm.true=NULL --> H does not exist

# If H exists and found something
if(!is.null(dgm$fs.harm.true) & !is.null(sg.harm.fs)){
any.H<-1.0
dfout<-df.pred
# PPV for H membership: Hpred --> treat.reccomend=0 & flag.harm=1
# Pred H and true H
aa<-with(dfout,sum(treat.recommend==0 & flag.harm==1))
# Pred Hc and true H
bb<-with(dfout,sum(treat.recommend==1 & flag.harm==1))
# Pred H and true Hc         
cc<-with(dfout,sum(treat.recommend==0 & flag.harm==0))
# Pred Hc and true Hc         
dd<-with(dfout,sum(treat.recommend==1 & flag.harm==0))

size.H<-with(dfout,sum(treat.recommend==0))
size.Hc<-with(dfout,sum(treat.recommend==1))

# hr estimate with true H
fit<-summary(coxph(cox.formula.sim,data=subset(df,flag.harm==1),robust=FALSE))$conf.int
hr.H.true<-c(fit[1])
l.H.true<-c(fit[3])
u.H.true<-c(fit[4])
rm("fit")

fit<-summary(coxph(cox.formula.sim,data=subset(df,flag.harm==0),robust=FALSE))$conf.int
hr.Hc.true<-c(fit[1])
l.Hc.true<-c(fit[3])
u.Hc.true<-c(fit[4])
rm("fit")

# hr estimate with estimated H
fit<-summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==0),robust=FALSE))$conf.int
hr.H.hat<-c(fit[1])
l.H.hat<-c(fit[3])
u.H.hat<-c(fit[4])
rm("fit")

fit<-summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==1),robust=FALSE))$conf.int
hr.Hc.hat<-c(fit[1])
l.Hc.hat<-c(fit[3])
u.Hc.hat<-c(fit[4])
rm("fit")

# Bias1 w.r.t hr.H.true
b1.H<-c(hr.H.hat-hr.H.true)
b2.H<-c(hr.H.hat-dgm$hr.H.true)

b1.Hc<-c(hr.Hc.hat-hr.Hc.true)
b2.Hc<-c(hr.Hc.hat-dgm$hr.Hc.true)

ppv<-aa/(aa+bb)
npv<-dd/(cc+dd)
specificity<-dd/(bb+dd)
sensitivity<-aa/(aa+cc)

#cat("PPV and NPV=",c(ppv,npv),"\n")
#cat("Sensitivity and Specificity=",c(specificity,sensitivity),"\n")

# Identify truth in dgm
found.1<-sum(grepl(dgm$fs.harm.true[1],sg.harm.fs))
found.2<-sum(grepl(dgm$fs.harm.true[2],sg.harm.fs))
found.both<-as.numeric(c(found.1+found.2)==2)
#found.al2<-as.numeric(length(sg.harm.fs)>=2)
found.al3<-as.numeric(length(sg.harm.fs)>=3)

}

# If H exists and did NOT find something
if(!is.null(dgm$fs.harm.true) & is.null(sg.harm.fs)){
  any.H<-0
  size.H<-0
  size.Hc<-nrow(df)
  
  # hr estimate with true H
#  hr.H.true<-c(summary(coxph(cox.formula.sim,data=subset(df,flag.harm==1),robust=FALSE))$conf.int[1])
#  hr.Hc.true<-c(summary(coxph(cox.formula.sim,data=subset(df,flag.harm==0),robust=FALSE))$conf.int[1])
  
  fit<-summary(coxph(cox.formula.sim,data=subset(df,flag.harm==1),robust=FALSE))$conf.int
  hr.H.true<-c(fit[1])
  l.H.true<-c(fit[3])
  u.H.true<-c(fit[4])
  rm("fit")
  
  fit<-summary(coxph(cox.formula.sim,data=subset(df,flag.harm==0),robust=FALSE))$conf.int
  hr.Hc.true<-c(fit[1])
  l.Hc.true<-c(fit[3])
  u.Hc.true<-c(fit[4])
  rm("fit")
  
 
  # hr estimate with estimated H, in this case, does not exist
  hr.H.hat<-NA
  l.H.hat<-NA
  u.H.hat<-NA
  
  # Define HC as ITT 
 # hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int[1])

 fit<-summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
  hr.Hc.hat<-c(fit[1])
  l.Hc.hat<-c(fit[3])
  u.Hc.hat<-c(fit[4])
  rm("fit")
  
   # Bias1 w.r.t hr.H.true
  b1.H<-NA
  b2.H<-NA
  
  b1.Hc<-c(hr.Hc.hat-hr.Hc.true)
  b2.Hc<-c(hr.Hc.hat-dgm$hr.Hc.true)
  
  # no H found, treat.recommend==1 for all
  # Pred H and true H
  #aa<-with(dfout,sum(treat.recommend==0 & flag.harm==1))
  # NO H found
  aa<-0
  # Pred Hc and true H
  #bb<-with(dfout,sum(treat.recommend==1 & flag.harm==1))
  # All treat.recommend=1
  bb<-with(df,sum(flag.harm==1))
  # Pred H and true Hc         
  #cc<-with(dfout,sum(treat.recommend==0 & flag.harm==0))
  cc<-0
  # Pred Hc and true Hc         
  #dd<-with(dfout,sum(treat.recommend==1 & flag.harm==0))
  dd<-with(df,sum(flag.harm==0))
  
  ppv<-aa/(aa+bb)
  npv<-dd/(cc+dd)
  specificity<-dd/(bb+dd)
  sensitivity<-ifelse(aa==0,0,aa/(aa+cc))
  
  # Identify truth in dgm
  found.1<-0
  found.2<-0
  found.both<-0
  found.al3<-0
  
}


# If H does NOT exist and found something
if(is.null(dgm$fs.harm.true) & !is.null(sg.harm.fs)){
  any.H<-1.0
  dfout<-df.pred
  # PPV for H membership: Hpred --> treat.reccomend=0 & flag.harm=1
  # Pred H and true H
  #aa<-with(dfout,sum(treat.recommend==0 & flag.harm==1))
  # flag.harm=0 for all
  aa<-0
  # Pred Hc and true H
  #bb<-with(dfout,sum(treat.recommend==1 & flag.harm==1))
  bb<-0
  # Pred H and true Hc         
  #cc<-with(dfout,sum(treat.recommend==0 & flag.harm==0))
  cc<-with(dfout,sum(treat.recommend==0))
  # Pred Hc and true Hc         
  #dd<-with(dfout,sum(treat.recommend==1 & flag.harm==0))
  dd<-with(dfout,sum(treat.recommend==1))
  
  size.H<-with(dfout,sum(treat.recommend==0))
  size.Hc<-with(dfout,sum(treat.recommend==1))
  
  # hr estimate with true H
  hr.H.true<-NA
  l.H.true<-NA
  u.H.true<-NA
  
  # no H --> ITT
  fit<-summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
  hr.Hc.true<-c(fit[1])
  l.Hc.true<-c(fit[3])
  u.Hc.true<-c(fit[4])
  rm("fit")
  
  
  # hr estimate with estimated H
  fit<-summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==0),robust=FALSE))$conf.int
  hr.H.hat<-c(fit[1])
  l.H.hat<-c(fit[3])
  u.H.hat<-c(fit[4])
  rm("fit")
  
  fit<-summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==1),robust=FALSE))$conf.int
  hr.Hc.hat<-c(fit[1])
  l.Hc.hat<-c(fit[3])
  u.Hc.hat<-c(fit[4])
  rm("fit")
  
  # Bias1 w.r.t hr.Hc.true
  b1.H<-NA
  b2.H<-NA
  
  b1.Hc<-c(hr.Hc.hat-hr.Hc.true)
  b2.Hc<-c(hr.Hc.hat-dgm$hr.Hc.true)
  
  ppv<-NA
  npv<-dd/(cc+dd)
  specificity<-dd/(bb+dd)
  sensitivity<-aa/(aa+cc)
  
  # Identify truth in dgm
  # fs.harm.true will be null
  # Still want to track how often v1.1 and v4.1 are found
  found.1<-sum(grepl(c("z1"),sg.harm.fs))
  found.2<-sum(grepl(c("z3"),sg.harm.fs))
  found.both<-as.numeric(c(found.1+found.2)==2)
  found.al3<-as.numeric(length(sg.harm.fs)>=3)
  
}

# If H does NOT exist and did NOT find something
# Define hat(Hc)=ITT
if(is.null(dgm$fs.harm.true) & is.null(sg.harm.fs)){
  any.H<-0
  # **flag.harm=0 for all**
  # Pred H and true H
  #aa<-with(dfout,sum(treat.recommend==0 & flag.harm==1))
  # flag.harm=0 for all
  aa<-0
  # Pred Hc and true H
  #bb<-with(dfout,sum(treat.recommend==1 & flag.harm==1))
  bb<-0
  # Pred H and true Hc         
  #cc<-with(dfout,sum(treat.recommend==0 & flag.harm==0))
  cc<-0
  # Pred Hc and true Hc         
  #dd<-with(dfout,sum(treat.recommend==1 & flag.harm==0))
  # Did NOT find anything --> treat.recommend=1 for all
  dd<-nrow(df)
  
  size.H<-0
  size.Hc<-nrow(df)
  
  # hr estimate with true H
  hr.H.true<-NA
  l.H.true<-NA
  u.H.true<-NA
  
  # ITT
  fit<-summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
  hr.Hc.true<-c(fit[1])
  l.Hc.true<-c(fit[3])
  u.Hc.true<-c(fit[4])
  rm("fit")
  
  
  # hr estimate with estimated H
  hr.H.hat<-NA
  l.H.hat<-NA
  u.H.hat<-NA
  
  #hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int[1])
  fit<-summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
  hr.Hc.hat<-c(fit[1])
  l.Hc.hat<-c(fit[3])
  u.Hc.hat<-c(fit[4])
  rm("fit")
  
  # Bias1 w.r.t hr.H.true
  b1.H<-NA
  b2.H<-NA
  
  b1.Hc<-c(hr.Hc.hat-hr.Hc.true)
  b2.Hc<-c(hr.Hc.hat-dgm$hr.Hc.true)
  
  ppv<-NA
  npv<-dd/(cc+dd)
  specificity<-dd/(bb+dd)
  sensitivity<-NA
  
  # Identify truth in dgm
  # fs.harm.true will be null
  # Still want to track how often v1.1 and v4.1 are found
  found.1<-sum(grepl(c("z1"),sg.harm.fs))
  found.2<-sum(grepl(c("z3"),sg.harm.fs))
  found.both<-as.numeric(c(found.1+found.2)==2)
  found.al3<-as.numeric(length(sg.harm.fs)>=3)
  
}

df.res<-data.table(any.H,size.H,size.Hc,ppv,npv,specificity,sensitivity,found.1,found.2,found.both,found.al3,
hr.H.true,hr.Hc.true,hr.H.hat,hr.Hc.hat,b1.H,b2.H,b1.Hc,b2.Hc,p.cens,analysis,taumax,
hr.itt,l.itt,u.itt,
hr.adj.itt,l.adj.itt,u.adj.itt,
l.H.true,u.H.true,
l.Hc.true,u.Hc.true,
l.H.hat,u.H.hat,
l.Hc.hat,u.Hc.hat)
return(df.res)
}







