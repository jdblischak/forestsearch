
count.weighted<-function(x,y,w=rep(1,length(y))){
  sum(w*(y<=x))
}

risk.weighted<-function(x,y,w=rep(1,length(y))){
  sum(w*(y>=x))
}

# ybar is risk process corresponding to timepoints 
# the default is at all (sorted) observed survival times 
# i.e., KM estimates at tpoints w.r.t. ybar(tpoints)
KM_estimates <- function(ybar,nbar){
  dN<-diff(c(0,nbar))
  dN.risk<-ifelse(ybar>0,dN/ybar,0.0)
  #chf <- cumsum(dN.risk)
  #var.chf<-cumsum(ifelse(ybar>0,dN/(ybar^2),0.0))
  S.KM <- cumprod(1-dN.risk)
  # return greenwood variance version
  aa<-dN
  bb<-ybar*(ybar-dN)
  var.KM<-(S.KM^2)*cumsum(ifelse(ybar>0,aa/bb,0.0))
  return(list(S.KM=S.KM,sig2.KM=var.KM))
}

wlr_estimates <- function(ybar0,ybar1,nbar0,nbar1,rho=0,gamma=0){
  # changing notation here
  dN.z0<- diff(c(0, nbar0))
  dN.z1<- diff(c(0, nbar1))
  dN.pooled<-dN.z0+dN.z1
  risk.z1 <- ybar1
  risk.z0 <- ybar0
  risk.pooled<-risk.z0+risk.z1
  ##########################
  # Pooled K-M estimator
  dN.Risk<-ifelse(risk.pooled>0,dN.pooled/risk.pooled,0)
  S.pool<-cumprod(1-dN.Risk)
  temp<-S.pool[1:(length(S.pool)-1)]
  S.pool<-c(1,temp) # S.pool(t-)
  w<-(S.pool^rho)*((1-S.pool)^gamma)
  K<-ifelse(risk.pooled>0,w*(risk.z0*risk.z1)/risk.pooled,0.0)
  term0<-sum(ifelse(risk.z0>0,(K/risk.z0)*dN.z0,0.0))
  term1<-sum(ifelse(risk.z1>0,(K/risk.z1)*dN.z1,0.0))
  lr<-term1-term0
  # variance
  h0<-ifelse(risk.z0==0,0,(K^2/risk.z0))
  h1<-ifelse(risk.z1==0,0,(K^2/risk.z1))
  dJ<-ifelse(risk.pooled==1,0,(dN.pooled-1)/(risk.pooled-1))
  dL<-ifelse(risk.pooled==0,0,dN.pooled/risk.pooled)
  sig2<-sum((h0+h1)*(1-dJ)*dL)
  return(list(lr=lr,sig2=sig2))
}

wlr_dhat_estimates <- function(dfcounting,rho=0,gamma=0,tzero=24){
  # changing notation here
  at.points <- dfcounting$at.points
  nbar0 <- dfcounting$nbar0
  nbar1 <- dfcounting$nbar1
  ybar0 <- dfcounting$ybar0
  ybar1 <- dfcounting$ybar1
  # KM estimates are already here
  S1 <- dfcounting$surv1 
  S0 <- dfcounting$surv0 
  Spool <- dfcounting$survP 
  # Estimates at tzero
  loc_tzero <- which.max(at.points > tzero)
  if(at.points[loc_tzero] <= tzero & at.points[loc_tzero+1] > tzero) 
  {
    dhat_tzero <- c(S1[loc_tzero]-S0[loc_tzero])
  }else{
    dhat_tzero <- c(S1[loc_tzero-1]-S0[loc_tzero-1])
  }
  # Common survival at tzero: S(tzero)
  Sp_tzero <- Spool[loc_tzero]
  dN.z0<- diff(c(0, nbar0))
  dN.z1<- diff(c(0, nbar1))
  dN.pooled<-dN.z0+dN.z1
  risk.z1 <- ybar1
  risk.z0 <- ybar0
  risk.pooled<-risk.z0+risk.z1
  ##########################
  # Pooled K-M estimator
  S.pool <- Spool 
  temp<-S.pool[1:(length(S.pool)-1)]
  S.pool<-c(1,temp) # S.pool(t-)
  # t- as weights need to be predictable (fixed prior to t)
  w<-(S.pool^rho)*((1-S.pool)^gamma)
  K<-ifelse(risk.pooled>0,w*(risk.z0*risk.z1)/risk.pooled,0.0)
  term0<-sum(ifelse(risk.z0>0,(K/risk.z0)*dN.z0,0.0))
  term1<-sum(ifelse(risk.z1>0,(K/risk.z1)*dN.z1,0.0))
  lr<-term1-term0
  # variance
  h0<-ifelse(risk.z0==0,0,(K^2/risk.z0))
  h1<-ifelse(risk.z1==0,0,(K^2/risk.z1))
  dJ<-ifelse(risk.pooled==1,0,(dN.pooled-1)/(risk.pooled-1))
  dL<-ifelse(risk.pooled==0,0,dN.pooled/risk.pooled)
  
  sig2_lr<-sum((h0+h1)*(1-dJ)*dL)
  
  # Term for covariance
  # w_integral_t0
  w_tzero <- w*c(ifelse(at.points <= tzero, 1,0))
  w_integral_t0 <- sum(w_tzero*(1-dJ)*dL)
  
  cov_wlr_dhat <- Sp_tzero*w_integral_t0
  
  # Term for variance of t0-difference 
  # assuming null
  
  h2 <- ifelse(risk.z0*risk.z1>0,(risk.pooled/(risk.z0*risk.z1)),0)
  h2 <- c(h2*ifelse(at.points <= tzero,1,0))
  
  sig2_dhat <- c(Sp_tzero^2)*sum(h2*(1-dJ)*dL)
  
  cor_wlr_dhat <- cov_wlr_dhat/c(sqrt(sig2_lr)*sqrt(sig2_dhat))
  return(list(lr=lr,sig2_lr=sig2_lr, dhat=dhat_tzero, cov_wlr_dhat=cov_wlr_dhat, sig2_dhat=sig2_dhat, cor_wlr_dhat=cor_wlr_dhat))
}

# data frame for counting process components
# returns risk, counting, KM, var(KM), and standard log-rank 
# note that these quantities are calculated at the (n=n0+n1) 
# sorted observed survival times
# If there is no strata in the dataset then set strata.name=NULL (default)
# which will turn-off stratified estimates

df_counting <- function(df,outcome.name,event.name,treat.name,strata.name=NULL){
  
  time <- df[,outcome.name]
  delta <- df[,event.name] 
  z <- df[,treat.name]  
  
  if(!is.null(strata.name))  strata <- df[,strata.name]
  if(is.null(strata.name)) strata <- rep("All", length(time))
  
  ans <- list()
  if(!all(z %in% c(0,1))) stop("Treatment must be numerical indicator: 0=control, 1=experimental")
  is.sorted<-!is.unsorted(time)
  if(!is.sorted){
    id <- order(time) 
    time <- time[id] 
    delta <- delta[id] 
    z<-z[id] 
    strata<- strata[id]
    stratum <- unique(strata)
  }
  # Calculate counting and at-risk process for overall and stratum-specific
  # at these timepoints
  # For stratum, the timepoints will be marked by their stratum membership
  # which for calculating stratified statistics, these would represent the observed
  # timepoints in the strata 
  # Note, here we use all time points -- not just unique
  # so that stratified estimates can be calculated where tied observations
  # can belong to different stratum
  at.points <- time
  # ybar0, nbar0 (overall)
  U0<-time[which(z==0)]
  D0<-delta[which(z==0)]
  # Risk and Counting processes
  ybar0 <-unlist(lapply(as.list(at.points),risk.weighted,y=U0))
  nbar0 <-unlist(lapply(as.list(at.points),count.weighted,y=U0[D0==1]))
  # KM and Var(KM) for z=0
  temp <- KM_estimates(ybar=ybar0,nbar=nbar0)
  surv0 <- temp$S.KM
  sig2_surv0 <- temp$sig2.KM
  # ybar1, nbar1 (overall)
  U1<-time[which(z==1)]
  D1<-delta[which(z==1)]
  # Risk and Counting processes
  ybar1 <-unlist(lapply(as.list(at.points),risk.weighted,y=U1))
  nbar1 <-unlist(lapply(as.list(at.points),count.weighted,y=U1[D1==1]))
  # KM and Var(KM) for z=0
  temp <- KM_estimates(ybar=ybar1,nbar=nbar1)
  surv1 <- temp$S.KM
  sig2_surv1 <- temp$sig2.KM
  # Pooled
  temp <- KM_estimates(ybar=c(ybar0+ybar1),nbar=c(nbar0+nbar1))
  survP <- temp$S.KM
  sig2_survP <- temp$sig2.KM
  # Log-rank 
  get_lr <- wlr_estimates(ybar0=ybar0,ybar1=ybar1,nbar0=nbar0,nbar1=nbar1,rho=0,gamma=0)
  ans$lr <- get_lr$lr
  ans$sig2_lr <- get_lr$sig2 
  
  ans$at.points <- at.points
  ans$strata <- strata
  # Control elements
  ans$ybar0 <- ybar0
  ans$nbar0 <- nbar0
  ans$surv0 <- surv0
  ans$sig2_surv0 <- sig2_surv0
  
  # Experimental elements
  ans$ybar1 <- ybar1
  ans$nbar1 <- nbar1
  ans$surv1 <- surv1
  ans$sig2_surv1 <- sig2_surv1
  
  ans$survP <- survP
  ans$sig2_survP <- sig2_survP
  
  # If stratified (strata >=2)
  # Calculate for each stratum
  # but return the processes corresponding to the strata value
  # That is if observation corresponding to at.points[1] had stratum=1
  # then return ybar0_strata = ybar[stratum=1]
  
  if(length(stratum)>1){
    ybar0_mat <- matrix(NA,nrow=length(at.points),ncol=length(stratum))
    nbar0_mat <- matrix(NA,nrow=length(at.points),ncol=length(stratum))
    ybar1_mat <- matrix(NA,nrow=length(at.points),ncol=length(stratum))
    nbar1_mat <- matrix(NA,nrow=length(at.points),ncol=length(stratum))
    
    surv0_mat <- matrix(NA,nrow=length(at.points),ncol=length(stratum))
    surv1_mat <- matrix(NA,nrow=length(at.points),ncol=length(stratum))
    sig2_surv0_mat <- matrix(NA,nrow=length(at.points),ncol=length(stratum))
    sig2_surv1_mat <- matrix(NA,nrow=length(at.points),ncol=length(stratum))
    
    survP_mat <- matrix(NA,nrow=length(at.points),ncol=length(stratum))
    sig2_survP_mat <- matrix(NA,nrow=length(at.points),ncol=length(stratum))
    
    lr_stratified <- 0.0
    sig2_lr_stratified <- 0.0
    
    for(ss in 1:length(stratum)){
      this_stratum <- stratum[ss]  
      U0_s<-time[which(z==0 & strata==this_stratum)]
      D0_s<-delta[which(z==0 & strata==this_stratum)]
      U1_s<-time[which(z==1 & strata==this_stratum)]
      D1_s<-delta[which(z==1 & strata==this_stratum)]
      # Risk and Counting processes
      
      ybar0_s <- unlist(lapply(as.list(at.points),risk.weighted,y=U0_s))
      nbar0_s <- unlist(lapply(as.list(at.points),count.weighted,y=U0_s[D0_s==1]))
      
      ybar1_s <- unlist(lapply(as.list(at.points),risk.weighted,y=U1_s))
      nbar1_s <- unlist(lapply(as.list(at.points),count.weighted,y=U1_s[D1_s==1]))
      
      ybar0_mat[, ss] <- ybar0_s 
      nbar0_mat[, ss] <- nbar0_s
      ybar1_mat[, ss] <- ybar1_s
      nbar1_mat[, ss] <- nbar1_s
      
      temp <- KM_estimates(ybar=ybar0_s,nbar=nbar0_s)
      surv0_mat[, ss] <- temp$S.KM
      sig2_surv0_mat[, ss] <- temp$sig2.KM
      rm("temp")
      temp <- KM_estimates(ybar=ybar1_s,nbar=nbar1_s)
      surv1_mat[, ss] <- temp$S.KM
      sig2_surv1_mat[, ss] <- temp$sig2.KM
      rm("temp")
      temp <- KM_estimates(ybar=c(ybar0_s+ybar1_s),nbar=c(nbar0_s+nbar1_s))
      survP_mat[, ss] <- temp$S.KM
      sig2_survP_mat[, ss] <- temp$sig2.KM
      rm("temp")
      temp <- wlr_estimates(ybar0=ybar0_s,ybar1=ybar1_s,nbar0=nbar0_s,nbar1=nbar1_s,rho=0,gamma=0)
      lr_stratified <- lr_stratified+temp$lr
      sig2_lr_stratified <- sig2_lr_stratified+temp$sig2
      rm("temp")
    }
    
    # For each strata value, extract the relevant column of ybar0_mat
    id_stratum <- unlist(lapply(strata,function(x){which(stratum==x)}))
    # Extract column corresponding to id_stratum
    index_extract <- cbind(1:length(at.points),id_stratum)
    
    ybar0_s <- ybar0_mat[index_extract]
    ybar1_s <- ybar1_mat[index_extract]
    nbar0_s <- nbar0_mat[index_extract]
    nbar1_s <- nbar1_mat[index_extract]
    
    surv0_s <- surv0_mat[index_extract]
    sig2_surv0_s <- sig2_surv0_mat[index_extract]
    
    surv1_s <- surv1_mat[index_extract]
    sig2_surv1_s <- sig2_surv1_mat[index_extract]
    
    survP_s <- survP_mat[index_extract]
    sig2_survP_s <- sig2_survP_mat[index_extract]
    
    ans$ybar0_stratum <- ybar0_s
    ans$nbar0_stratum <- nbar0_s
    ans$surv0_stratum <- surv0_s
    ans$sig2_surv0_stratum <- sig2_surv0_s
    
    ans$ybar1_stratum <- ybar1_s
    ans$nbar1_stratum <- nbar1_s
    ans$surv1_stratum <- surv1_s
    ans$sig2_surv1_stratum <- sig2_surv1_s
    
    ans$survP_stratum <- survP_s 
    ans$sig2_survP_stratum <- sig2_survP_s
    ans$lr_stratified <- lr_stratified
    ans$sig2_lr_stratified <- sig2_lr_stratified 
  }
  return(ans)
}

# For given at.points we want KM estimates for both arms at the same timepoints (at.points)
KM_diff <- function(df,outcome.name,event.name,treat.name,at.points=c(sort(df[,outcome.name])),alpha=0.05){
  time <- df[,outcome.name]
  delta <- df[,event.name] 
  z <- df[,treat.name]  
  ans <- list()
  if(!all(z %in% c(0,1))) stop("Treatment must be numerical indicator: 0=control, 1=experimental")
  is.sorted<-!is.unsorted(time)
  if(!is.sorted){
    id <- order(time) 
    time <- time[id] 
    delta <- delta[id] 
    z<-z[id] 
  }
  # Calculate counting and at-risk process for overall and stratum-specific
  ans$at.points <- at.points
  U0<-time[which(z==0)]
  D0<-delta[which(z==0)]
  # Risk and Counting processes
  ybar0 <-unlist(lapply(as.list(at.points),risk.weighted,y=U0))
  nbar0 <-unlist(lapply(as.list(at.points),count.weighted,y=U0[D0==1]))
  # KM and Var(KM) for z=0
  temp <- KM_estimates(ybar=ybar0,nbar=nbar0)
  surv0 <- temp$S.KM
  sig2_surv0 <- temp$sig2.KM
  ans$surv0 <- surv0
  ans$sig2_surv0 <- sig2_surv0
  U1<-time[which(z==1)]
  D1<-delta[which(z==1)]
  # Risk and Counting processes
  ybar1 <-unlist(lapply(as.list(at.points),risk.weighted,y=U1))
  nbar1 <-unlist(lapply(as.list(at.points),count.weighted,y=U1[D1==1]))
  # KM and Var(KM) for z=0
  temp <- KM_estimates(ybar=ybar1,nbar=nbar1)
  surv1 <- temp$S.KM
  sig2_surv1 <- temp$sig2.KM
  ans$surv1 <- surv1
  ans$sig2_surv1 <- sig2_surv1
  # Difference
  dhat <- surv1-surv0
  sig2_dhat <- sig2_surv0+sig2_surv1
  ans$dhat <- dhat
  ans$sig2_dhat <- sig2_dhat
  # lower and upper bounds (point-wise)
  c_alpha <- qnorm(1-alpha/2)
  lower <- dhat - c_alpha*sqrt(sig2_dhat)
  upper <- dhat + c_alpha*sqrt(sig2_dhat)
  ans$lower <- lower
  ans$upper <- upper
  return(ans)
}


# For subgroup displays which includes ITT
plotKM.band_subgroups<-function(df,outcome.name,event.name,treat.name,
                                sg_labels=c("cps_10","cps_5"),ltype="s",lty=1,draws=20,lwd=2,
                                sg_colors=c("blue","brown"),ymax.pad=0.0,ymin.pad=0.0,
                                taus=c(-Inf,Inf),yseq_length=5,cex_Yaxis=0.8,risk_cex=0.8,
                                by.risk=6,risk.add=NULL,xmax=NULL,ymin=NULL,ymax=NULL, ymin.del=0.035,
                                y.risk1=NULL, y.risk2=NULL,ymin2=NULL, risk_offset=NULL, risk.pad=0.01,
                                risk_delta=0.0275,tau_add=NULL,time.zero.pad=0, time.zero.label=0.0,
                                xlabel=NULL,ylabel=NULL,color="lightgrey",Maxtau=NULL,
                                ylim=c(min(lower,na.rm=TRUE),max(upper,na.rm=TRUE))){
  
  if(length(sg_labels) != length(sg_colors)) stop("SG labels and colors do not match")  
  
  if(is.null(risk_offset)) risk_offset = (1+length(sg_labels))*risk_delta
  
  time.zero <- time.zero.label-time.zero.pad
  
  if(is.null(ylabel)) ylabel=c(expression(delta(t)==hat(S)[1](t)-hat(S)[0](t)))
  if(is.null(xlabel)) xlabel="Months"  
  
  # ITT 
  Y<-df[,outcome.name]
  E<-df[,event.name]
  Treat<-df[,treat.name]
  
  maxe_0 <- max(Y[E==1 & Treat==0])
  maxe_1 <- max(Y[E==1 & Treat==0])
  max_tau <- min(c(maxe_0,maxe_1))
  if(is.null(Maxtau)) max_tau <- max(c(max_tau,tau_add))
  if(!is.null(Maxtau)) max_tau <- Maxtau
  
  # add by.risk points
  riskpoints <- seq(0,max_tau,by=by.risk)
  
  # Truncate the observed timepoints at max_tau
  atpoints <- c(Y[which(Y<=max_tau)])
  at.points <- sort(unique(c(atpoints,max_tau,riskpoints)))
  
  fit <- KM_diff(df=df,outcome.name=outcome.name,event.name=event.name,
  treat.name=treat.name,at.points=at.points,alpha=0.05)
  
  xx0 <- fit$at.points
  yy0 <- fit$dhat
  l0 <- fit$lower
  u0 <- fit$upper
  
  at.points <- xx0
  risk.points<-round(c(seq(time.zero.label,max(at.points),by=by.risk)))
  risk.points<-sort(unique(c(risk.points,risk.add)))
  
  # Here time.zero artificially adds time buffer
  risk.points <- c(time.zero,risk.points) 
  
  # Total risk for ITT
  risk0 <-unlist(lapply(as.list(risk.points),risk.weighted,y=Y))
  
  risk0<-round(risk0)
  
  # Subgroups (SG) defined by sg_labels
  # Note: SG subgroup KM differences are calculated at same time points as ITT
  # The same time points are taus=(0,max_tau) (for all subgroups)
  # So store in matrix for plotting
  
  # Store KM differences for subgroups in Dsg_mat
  Dsg_mat <- matrix(NA,nrow=length(xx0),ncol=length(sg_labels))
  # Risk sets in Rsg_mat
  Rsg_mat <- matrix(NA,nrow=length(risk0),ncol=length(sg_labels))
  
  # Store KM estimates for control S0_mat
  S0_mat <- matrix(NA,nrow=length(xx0),ncol=length(sg_labels))
  S1_mat <- matrix(NA,nrow=length(xx0),ncol=length(sg_labels))
  
  for(sg in 1:length(sg_labels)){
    sg_flag <- sg_labels[sg]
    # In case missing sg_flag
    dfs <- df[,c(outcome.name,event.name,treat.name,sg_flag)]
    dfs <- na.omit(dfs)
    df_sg<-subset(dfs,dfs[,sg_flag]==1)
    
    Y_sg<-df_sg[,outcome.name]
    E_sg<-df_sg[,event.name]
    Treat_sg<-df_sg[,treat.name]
    
    res <- KM_diff(df=df_sg,outcome.name=outcome.name,event.name=event.name,
                   treat.name=treat.name,at.points=at.points,alpha=0.05)
    
    
    Dsg_mat[,sg] <- res$dhat
    
    rr <-unlist(lapply(as.list(risk.points),risk.weighted,y=Y_sg))
    
    Rsg_mat[,sg] <- c(round(rr))
    
    S0_mat[,sg] <- res$surv0
    S1_mat[,sg] <- res$surv1
  }
  
  x <- xx0
  mean.value <- yy0
  lower <- l0
  upper <- u0
  
  if(is.null(ymax)){
  ymax <- max(Dsg_mat,na.rm=TRUE)+ymax.pad
  ymax <- round(max(c(u0,ymax),na.rm=TRUE),2)
  }
  
  if(is.null(ymin)){
    ymin <- min(Dsg_mat,na.rm=TRUE)-ymin.pad
    ymin <- round(min(c(l0,ymin),na.rm=TRUE),2)
  }
  
  if(is.null(ymin2)) ymin2 <- ymin-risk_offset  

  yrisks <- seq(ymin2,ymin-(ymin.del+risk.pad),length=1+length(sg_labels))
  yrisks <- yrisks[order(yrisks,decreasing=TRUE)]
  
  plot(x[order(x)],mean.value[order(x)],type="n",axes=FALSE,xlab=xlabel,lty=lty,
  ylab=ylabel,ylim=c(ymin2,ymax),cex.lab=cex_Yaxis)
  polygon(c(x[order(x)],rev(x[order(x)])),c(lower[order(x)],rev(upper[order(x)])),col=color,border=FALSE)
  lines(x[order(x)],mean.value[order(x)],lty=lty,lwd=lwd,type=ltype)
  
  # horizontal line just below ymin 
  abline(h=ymin-ymin.del,lty=1,col="black")
  abline(h=time.zero.label,lty=2,col="black",lwd=0.5)
  
  for(sg in 1:length(sg_labels)){
    lines(x,Dsg_mat[,sg],col=sg_colors[sg],type=ltype,lwd=lwd)
  }
  
  if((ymax-ymin2) <= 0.5) ypoints <- round(seq(ymin,ymax,length=yseq_length),2)
  if((ymax-ymin2) > 0.5) ypoints <- seq(ymin,ymax,by=0.15)
  
  ypoints <- sort(c(time.zero.label,ypoints))
  
  axis(2,at=c(ypoints),cex.axis=cex_Yaxis, las=1)
  
  risk.points.label <- as.character(c(time.zero.label,risk.points[-1]))
  
  axis(1,at=c(risk.points),labels=c(risk.points.label),cex.axis=cex_Yaxis)
  
  box()
  # First, ITT risk
  text(c(risk.points),c(yrisks[1]),c(risk0),col="black",cex=risk_cex)
  for(sg in 1:length(sg_labels)){
    text(c(risk.points),c(yrisks[sg+1]),c(Rsg_mat[,sg]),col=sg_colors[sg],cex=risk_cex)
  }
  
  return(list(xpoints=xx0,Dhat_subgroups=Dsg_mat,s0_subgroups=S0_mat,s1_subgroups=S1_mat,
              rpoints=risk.points,Risk_subgroups=Rsg_mat,mean=mean.value,lower=lower,upper=upper))  
  
}


# General forestplot
SG_forestplot <-  function(df,arrow_text=c("favors Treatment","Control"),E.name=c("Treat"),C.name=c("Control"),outcome.name,event.name,treat.name,
                           sg_labels=c("Age<=65","Age>65"),footnote_text=NULL){
    title_text <- NULL
    # ITT 
    res_itt <- SG_HRtable2(dfa=df,df.OOB=NULL,fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,
                           sg_name="ITT",E.name=E.name,C.name=C.name)
    # Add separator header 
    zz <- res_itt  
    zz$Subgroup <- "                 "
    zz[,c(E.name,C.name)] <- ""
    zz[,c("est","low","hi","se")] <- NA
}

