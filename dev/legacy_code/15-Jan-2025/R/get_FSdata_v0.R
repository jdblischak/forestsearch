
# Setting up dataset for ForestSearch
# which includes options to apply lasso for dimension reduction
# and GRF cuts

is.continuous<-function(x,cutoff=4){ifelse(length(unique(x))>=cutoff,1,2)}

qlow <- function(x) c(quantile(x,0.25))
qhigh <- function(x) c(quantile(x,0.75))

# For continuous variables in conf_force_names
# setup for mean, median, qlow, and qhigh

cut_var <- function(x){
mx <- paste0("mean(",x,")")  
a <- paste0(x," <= ",mx)
mdx <- paste0("median(",x,")")  
b <- paste0(x," <= ",mdx)
qlx <- paste0("qlow(",x,")")  
c <- paste0(x," <= ",qlx)
qhx <- paste0("qhigh(",x,")")  
d <- paste0(x," <= ",qhx)
return(c(a,b,c,d))
}

get_conf_force <- function(df,conf.force.names,cont.cutoff=10){
res <- NULL
for(ccs in 1:length(conf.force.names)){
aa<-df[,c(conf.force.names[ccs])]
flag_cont <- is.continuous(aa,cutoff=cont.cutoff)
if(flag_cont==1){
# create mean, median, qlow, and qhigh  
this_cut <- cut_var(x=conf.force.names[ccs])
res <- c(res,this_cut)
}
rm("aa")
}
return(res)
}

get_FSdata<-function(df,use_lasso=FALSE,use_grf=FALSE,use_grf_only=FALSE,grf_cuts,confounders.name,
cont.cutoff=4,conf_force=NULL, conf.cont_medians=NULL, conf.cont_medians_force=NULL,
replace_med_grf=TRUE, defaultcut_names=NULL, cut_type="default", exclude_cuts = NULL,
outcome.name=NULL, event.name=NULL, details=TRUE){

# Default cuts forced per defaultcut_names
if(!is.null(defaultcut_names)){
conf_force_default <- get_conf_force(df=df,conf.force.names=defaultcut_names,cont.cutoff=4)
# append to conf_force
conf_force <- c(conf_force,conf_force_default)
}  
  if(use_grf & is.null(grf_cuts)){
    warning("GRF cuts need to be specified (run prior): Considering no cuts per GRF turning off GRF")
    use_grf<-FALSE
  }  
  if(use_lasso &  (is.null(outcome.name) | is.null(event.name))) stop("Cox variable names needed for lasso")
  
  # Flag continuous and categorical (continuous if unique values exceeds con.cutoff)
  flag_continuous<-rep(NA,length(confounders.name))
  for(ccs in 1:length(confounders.name)){
    aa<-df[,c(confounders.name[ccs])]
    flag_cont <- is.continuous(aa,cutoff=cont.cutoff)
    flag_continuous[ccs]<-(flag_cont==1)  
    rm("aa")
  }
  
if(details){
cat("# of continuous/categorical characteristics",c(sum(flag_continuous),sum(!flag_continuous)),"\n")
if(sum(flag_continuous)>0)  cat("Continuous characteristics:",c(confounders.name[flag_continuous]),"\n")
if(sum(!flag_continuous)>0)   cat("Categorical characteristics:",c(confounders.name[!flag_continuous]),"\n")
}
  if(sum(flag_continuous)==0) conf.categorical <- confounders.name
  if(sum(flag_continuous) > 0) conf.categorical <- confounders.name[!flag_continuous]
  
  if(is.null(conf.cont_medians) & is.null(conf.cont_medians_force) & sum(flag_continuous)>0){
  conf.cont_medians<-confounders.name[flag_continuous]
  }  
  
  lassokeep <- NULL
  lassoomit <- NULL
  
  if(use_lasso){
    set.seed(8316951)
    # Reduce dimension via Cox lasso
    xx<-as.matrix(df[,confounders.name])
    yy<-as.matrix(df[,c(outcome.name,event.name)])
    colnames(yy)<-c("time","status")
    
    cvfit <- cv.glmnet(xx, yy, family="cox") #first do 10-fold cross-validation to select lambda
    
  if(details) cat("CV lambda =",c(cvfit$lambda.min),"\n")
    
    m <- glmnet(xx, yy, family="cox", lambda=cvfit$lambda.min) #plugin the optimal lambda
    
    lassokeep<-confounders.name[which(m$beta!=0)]
    lassoomit<-confounders.name[which(m$beta==0)]
    
    if(details){
      print(m$beta)
      cat("Cox-LASSO selected:",c(lassokeep),"\n")
      cat("Cox-LASSO not selected:",c(lassoomit),"\n")
    }
    
    
if(length(lassokeep)>0){
   if(!is.null(conf.cont_medians)){
      to_keep<-which(conf.cont_medians %in% lassokeep)
      if(length(to_keep)>0) conf.cont_medians <- conf.cont_medians[to_keep]
      if(length(to_keep)==0) conf.cont_medians <- NULL
      }
    if(!is.null(conf.categorical)){
      to_keep<-which(conf.categorical %in% lassokeep)
      if(length(to_keep)>0)  conf.categorical <- conf.categorical[to_keep]
      if(length(to_keep)==0) conf.categorical <- NULL
      }
    } # If any selected per Lasso
  } # Done Lasso
 
if(use_grf_only) conf.cont_medians <- NULL
# If forcing any cuts, then done below
  
if(use_lasso & cut_type=="default"){
conf_force_lasso <- NULL
# Override conf.cont_medians
# conf.cont_medians are continuous factors selected per lasso  
# Create default cuts for conf.cont_medians not in defaultcut_names  
if(!is.null(defaultcut_names)){  
nottocut <- conf.cont_medians %in% defaultcut_names
lasso_tocut <- conf.cont_medians[!nottocut]    
if(length(lasso_tocut)>0) conf_force_lasso <- get_conf_force(df=df,conf.force.names=lasso_tocut,cont.cutoff=4)
# Override cuts at medians
conf.cont_medians <- NULL
}
if(is.null(defaultcut_names)){  
    lasso_tocut <- conf.cont_medians    
if(length(lasso_tocut)>0) conf_force_lasso <- get_conf_force(df=df,conf.force.names=lasso_tocut,cont.cutoff=4)
    # Override cuts at medians
    conf.cont_medians <- NULL
  }
# Append to conf_force
conf_force <- c(conf_force,conf_force_lasso)
if(details){
cat("Default cuts included from Lasso:",c(conf_force_lasso),"\n")
cat("Categorical after Lasso:",c(conf.categorical),"\n")
}  
 }

if(use_lasso & cut_type=="median"){
    if(!is.null(defaultcut_names)){  
      nottocut <- conf.cont_medians %in% defaultcut_names
      conf.cont_medians <- conf.cont_medians[!nottocut]    
    }
    if(details){
      cat("Median cuts included from Lasso:",c(conf.cont_medians),"\n")
      cat("Categorical after Lasso:",c(conf.categorical),"\n")
    }  
  }

  
  if(!use_lasso & cut_type=="default"){
    conf_force_add <- NULL
    # Override conf.cont_medians
    # conf.cont_medians are continuous factors selected per lasso  
    # Create default cuts for conf.cont_medians not in defaultcut_names  
    if(!is.null(defaultcut_names)){  
      nottocut <- conf.cont_medians %in% defaultcut_names
      tocut <- conf.cont_medians[!nottocut]    
if(length(tocut)>0) conf_force_add <- get_conf_force(df=df,conf.force.names=tocut,cont.cutoff=4)
      # Override cuts at medians
      conf.cont_medians <- NULL
    }
    if(is.null(defaultcut_names)){  
      tocut <- conf.cont_medians    
if(length(tocut)>0) conf_force_add <- get_conf_force(df=df,conf.force.names=tocut,cont.cutoff=4)
      # Override cuts at medians
      conf.cont_medians <- NULL
    }
    # Append to conf_force
    conf_force <- c(conf_force,conf_force_add)
    if(details){
      toprint <- min(20,length(conf_force_add))
      cat("Default cuts included (1st 20)",c(conf_force_add[1:toprint]),"\n")
      cat("Categorical:",c(conf.categorical),"\n")
    }  
  }
  
  if(!use_lasso & cut_type=="median"){
    if(!is.null(defaultcut_names)){  
      nottocut <- conf.cont_medians %in% defaultcut_names
      conf.cont_medians <- conf.cont_medians[!nottocut]    
    }
    if(details){
      toprint <- min(20,length(conf.cont_medians))
      cat("Median cuts included:",c(conf.cont_medians[1:toprint]),"\n")
      cat("Categorical:",c(conf.categorical),"\n")
    }  
  }

if(!is.null(conf.cont_medians_force)){
    conf.cont_medians<-c(conf.cont_medians,conf.cont_medians_force)
  }

if(details & use_grf){
cat("Factors per GRF:",c(grf_cuts),"\n")
}    

if(details & use_grf & length(conf.cont_medians)>0){
toprint <- min(20,length(conf.cont_medians))
cat("Continuous factors initially cut at medians:",c(conf.cont_medians[1:toprint]),"\n")
}    

# Remove any factors to cut at median if already in GRF
if(replace_med_grf){  
if(use_grf & length(conf.cont_medians)>0){
if(length(grf_cuts)>0){  
    # Which to omit  
    cMed_flag<-NULL
    for(ccs in 1:length(grf_cuts)){
    thiscut<-grf_cuts[ccs]
    flag_omit<-c(unlist(lapply(conf.cont_medians,function(x){grepl(x,thiscut)})))
  if(any(flag_omit))  cMed_flag<-c(cMed_flag,c(conf.cont_medians[flag_omit]))
   }  
if(length(cMed_flag)>0){
to_exclude <- (conf.cont_medians %in% cMed_flag)
# If excluding at least 1, the following remain
if(sum(!to_exclude)>0)  conf.cont_medians <- conf.cont_medians[!to_exclude]
# If excluding all then set to null
if(sum(!to_exclude)==0)  conf.cont_medians <- NULL
     } # Otherwise, NO factors omitted
   }
} 
if(details & use_grf & length(conf.cont_medians)>0)  cat("Factors after removing any duplicates also in GRF:",c(conf.cont_medians),"\n")
}

if(details & cut_type=="median" & use_lasso & length(conf.cont_medians)==0){
cat("***conf.cont_medians is NULL --> NO MEDIAN CUTS per lasso***","\n")
}  

# Re-introduce conf.cont_force_medians
if(!is.null(conf.cont_medians_force)) conf.cont_medians <- c(conf.cont_medians,conf.cont_medians_force)

  conf.cont_Medcuts<-NULL

    if(!is.null(conf.cont_medians)){
    # Default is to specify median cuts
    # Continuous covariate indexes
    for(aa in 1:length(conf.cont_medians)){
      mval <- as.character(round(median(df[,c(conf.cont_medians[aa])]),2))
      cc1<-paste0(conf.cont_medians[aa]," <= ",sep="")  
      cc2<-paste0(cc1,mval, sep="")
      conf.cont_Medcuts<-c(conf.cont_Medcuts,cc2)
       }
    }
  
  confs<-c(conf.categorical,conf.cont_Medcuts)

  # At this stage, these are confs per Lasso (GRF step is next)
  if(use_lasso) confs_lasso <- confs 
  if(use_grf){
  # If element in grf_cuts matches conf.categorical
  # then this is redundant per FS handling of categorical
  # factors (Any level will already be incorporated)
  # So remove grf_cuts already in conf.categorical
  if(details)  cat("Initial GRF cuts included",c(grf_cuts),"\n")  
  if(length(confs)>0 & length(grf_cuts)>0){  
  grf_cuts_keep<-NULL
  for(ccs in 1:length(grf_cuts)){
   thiscut<-grf_cuts[ccs]
  # Covariate corresponding to
  flag_omit<-c(unlist(lapply(confs,function(x){grepl(x,thiscut)})))
  # Only keep if does not already exist in conf.categorical 
  if(all(!flag_omit)) grf_cuts_keep<-c(grf_cuts_keep,thiscut)
  }  
  
  if(length(confs)==0 & length(grf_cuts)>0){  
  grf_cuts_keep <- grf_cuts
  }
  confs <- unique(c(confs,grf_cuts_keep))
  }
  } # use_grf
  # Factors included per GRF not in confs_lasso

if(use_lasso & use_grf){
which_both <- (confs %in% confs_lasso)
if(details & length(which_both) >0){
cat("Factors included per GRF (not in lasso)",c(confs[!which_both]),"\n")
}
  }
  
if(!is.null(conf_force)){
conf_forceNew <- NULL
for(cfs in 1:length(conf_force)){
loc_name <- charmatch(confounders.name,conf_force[cfs])
# return correct name (e.g, differentiate x1, x10)
index_name <- which(loc_name==1)
# Find exact match
if(length(index_name)>1){
  confs2 <- confounders.name[index_name]
  lc <- str_length(confs2)
  ctoget <- str_sub(conf_force[cfs],1,max(lc))
   # which to extract
  itoget <- which(confounders.name == ctoget)
   cfs_name <- confounders.name[itoget]
}
if(length(index_name)==1) cfs_name <- confounders.name[which(loc_name==1)]
## mean, median, qlow, or qhigh?
calcmean <- grepl("mean",conf_force[cfs])
calcmedian <- grepl("median",conf_force[cfs])
calcqlow <- grepl("qlow",conf_force[cfs])
calcqhigh <- grepl("qhigh",conf_force[cfs])
if(calcmean){
  mval <- as.character(round(mean(df[,cfs_name]),1))
  cc1<-paste0(cfs_name," <= ",sep="")  
  cc2<-paste0(cc1,mval, sep="")
  conf_forceNew<-c(conf_forceNew,cc2)
}
if(calcmedian){
  mval <- as.character(round(median(df[,cfs_name]),1))
  cc1<-paste0(cfs_name," <= ",sep="")  
  cc2<-paste0(cc1,mval, sep="")
  conf_forceNew<-c(conf_forceNew,cc2)
}
if(calcqlow){
  mval <- as.character(round(quantile(df[,cfs_name],0.25),1))
  cc1<-paste0(cfs_name," <= ",sep="")  
  cc2<-paste0(cc1,mval, sep="")
  conf_forceNew<-c(conf_forceNew,cc2)
}

if(calcqhigh){
  mval <- as.character(round(quantile(df[,cfs_name],0.75),1))
  cc1<-paste0(cfs_name," <= ",sep="")  
  cc2<-paste0(cc1,mval, sep="")
  conf_forceNew<-c(conf_forceNew,cc2)
}
if(!calcmean & !calcmedian & !calcqlow & !calcqhigh){
conf_forceNew<-c(conf_forceNew,conf_force[cfs])
}
  }
  }

  
if(!is.null(conf_force)) confs<-unique(c(confs,conf_forceNew)) 

# Excluding cuts
if(!is.null(exclude_cuts)){  
# Remove the restricted cuts (eg., not allowing a variable to be "<=")
for(ee in 1:length(exclude_cuts)){  
to_exclude <- grepl(exclude_cuts[ee],confs)
confs <- confs[!to_exclude]
}
 } 
  n_confs<-length(confs)
  if(n_confs==0) stop("Error in FS dataset prior to flag drop")
  
  # Arrange by continuous and categorical
  # Classify as continuous or categorical
  flag_continuous<-rep(NA,length(confs))
  flag_drop<-rep(NA,length(confs)) # Flag factors with only 1 level and drop
  for(ccs in 1:length(confs)){
    thiscut<-confs[ccs]
    # Covariate corresponding to
    cov_index<-c(unlist(lapply(confounders.name,function(x){grepl(x,thiscut)})))
    cut_name<-confounders.name[cov_index]
    aa<-df[,c(cut_name)]
    flag_continuous[ccs]<-(is.continuous(aa,cutoff=cont.cutoff)==1)  
    flag_drop[ccs]<-(length(unique(aa))<=1)
  }
  
  if(details & any(flag_drop)){
    cat("Dropping variables (only 1 level):",c(confs[flag_drop]),"\n")
  }
  
  conf.categorical<-confs[!flag_continuous & !flag_drop]
  
  conf.cont_cuts <- NULL
  if(sum(flag_continuous)>0){
  conf.cont_cuts<-confs[flag_continuous & !flag_drop]
  }
  
  # Note, possibility that a cut point only has 1 level
  # Drop any such factors

flag_drop2<-NULL    
if(length(conf.cont_cuts)>=1){
  # Note, need to copy (otherwise qcheck will not be saved) 
    df_check<-copy(df)
  flag_drop2<-rep(NA,length(conf.cont_cuts))
    for(ccs in 1:length(conf.cont_cuts)){
      thiscut<-conf.cont_cuts[ccs]
      df_check <- within(df_check,{
        qcheck <- as.factor(ifelse(eval(parse(text=thiscut)),1,0))
      })
      flag_drop2[ccs]<-(length(unique(df_check$qcheck))<=1)
    }
  
  if(details & any(flag_drop2)){
    cat("Dropping variables (cut only has 1 level)",c(conf.cont_cuts[flag_drop2]),"\n")
  }
conf.cont_cuts<-conf.cont_cuts[!flag_drop2]
rm("df_check")
  } # flag_drop2   
  
# n_confs after drop1 and drop2
n_confs<-n_confs-(sum(flag_drop)+sum(flag_drop2))

confs<-c(conf.cont_cuts,conf.categorical)

if(length(confs)!=n_confs) stop("Error in categorical and continuous classifications")
  
  # These will be names q1,..,qK (say)
  names_new<-c(unlist(lapply(c(1:length(confs)),function(x){paste0("q",x,sep="")})))
  
  # For the continuous
  if(length(conf.cont_cuts)>=1){
    for(ccs in 1:length(conf.cont_cuts)){
      thiscut<-conf.cont_cuts[ccs]
      df<-within(df,{
        q <- as.factor(as.numeric(ifelse(eval(parse(text=thiscut)),1,0)))
      })
      loc_names<-which(names(df)=="q")
      names(df)[loc_names]<-names_new[ccs]
    }
  }
  
  if(length(conf.categorical)>=1){
    for(ccs in 1:length(conf.categorical)){
      thiscut<-conf.categorical[ccs]
      df<-within(df,{
        q <- as.factor(as.numeric(eval(parse(text=thiscut))))
      })
      loc_names<-which(names(df)=="q")
      names(df)[loc_names]<-names_new[ccs+length(conf.cont_cuts)] # Add cont index 
    }
  }

# NOTE: include check that all confs are 0,1 (not TRUE, FALSE)
check_factors<-apply(df[,c(names_new)],2,function(x){c(unique(x))})
if(any(c(check_factors) %in% c("TRUE","FALSE"))) stop("Error in factor setup")
  
  if(details){
    cat("# of candidate subgroup factors=",c(length(confs)),"\n")
    print(confs)
  }  
  
  return(list(df=df,confs_names=names_new,confs=confs,lassokeep=lassokeep,lassoomit=lassoomit))
}


find_qloc <- function(x){
# First 9
loc9 <- substr(x,3,3)
loc10 <- substr(x,4,4)
loc100 <- substr(x,5,5)
loc1000 <- substr(x,6,6)
if(loc9==".") qloc <- "first9"
if(loc10==".") qloc <- "tenplus"
if(loc100==".") qloc <- "hundredplus"
if(loc1000==".") qloc <- "thousandplus"
return(qloc)
}

FS_labels<-function(Qsg,confs_labels){
  label_index<-c(1:length(confs_labels))
  nfacs <- length(Qsg)
  sg_labels <- NULL
  for(sgs in 1:nfacs){
  # Check how many elements
  # More than 4 --> factor is 10th or greater 
  qsg <- Qsg[sgs]
  # Find "." for location of action (e.g., q1.0)
  # Determine whether 
  # first 9 "first9" 
  # 10-99 "tenplus" or 
  # 100+ "hundredplus"
  qloc <- find_qloc(qsg)
  #flag_4<-as.numeric(substr(qsg,5,5))
  # If 5th element is missing then factor is 9th or less
  if(qloc=="first9"){
    # extract index from qsg
    qindx<-as.numeric(substr(qsg,2,2)) # index <= 9
    qaction<-as.numeric(substr(qsg,4,4))
    sg_label<-confs_labels[qindx]
    if(qaction==0) sg_label <- paste0("!{",sg_label,sep="}")
    if(qaction==1) sg_label <- paste0("{",sg_label,sep="}")
  }
  if(qloc=="tenplus"){
    # extract index from qsg
    qindx<-as.numeric(substr(qsg,2,3))
    qaction<-as.numeric(substr(qsg,5,5))
    sg_label<-confs_labels[qindx]
    if(qaction==0) sg_label <- paste0("!{",sg_label,sep="}")
    if(qaction==1) sg_label <- paste0("{",sg_label,sep="}")
  }
  if(qloc=="hundredplus"){
    # extract index from qsg
    qindx<-as.numeric(substr(qsg,2,4))
    qaction<-as.numeric(substr(qsg,6,6))
    sg_label<-confs_labels[qindx]
    if(qaction==0) sg_label <- paste0("!{",sg_label,sep="}")
    if(qaction==1) sg_label <- paste0("{",sg_label,sep="}")
  }
  if(qloc=="thousandplus"){
    # extract index from qsg
    qindx<-as.numeric(substr(qsg,2,5))
    qaction<-as.numeric(substr(qsg,7,7))
    sg_label<-confs_labels[qindx]
    if(qaction==0) sg_label <- paste0("!{",sg_label,sep="}")
    if(qaction==1) sg_label <- paste0("{",sg_label,sep="}")
  }
  sg_labels <- c(sg_labels,sg_label)
  }
  return(sg_labels)
}
