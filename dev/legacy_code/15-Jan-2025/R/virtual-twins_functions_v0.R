
formatRCTDataset2<-function (dataset, outcome.field, treatment.field, interactions = TRUE) 
{
  if (!is.data.frame(dataset)) 
    stop("Dataset parameter must be data.frame")
  if (!is.character(outcome.field)) 
    stop(sprintf("%s, outcome.field parameter must be a string", 
                 outcome.field))
  if (!is.character(treatment.field)) 
    stop(sprintf("%s, treatment.field parameter must be a string", 
                 treatment.field))
  if (!outcome.field %in% colnames(dataset)) 
    stop(sprintf("%s must be in data.frame colnames", outcome.field))
  outcome.field.which <- which(outcome.field == colnames(dataset))
  if (!treatment.field %in% colnames(dataset)) 
    stop(sprintf("%s must be in data.frame colnames", treatment.field))
  treatment.field.which <- which(treatment.field == colnames(dataset))
  d <- dataset
  outcome <- d[, outcome.field.which]
  if (!is.factor(outcome)) 
    outcome <- as.factor(outcome)
  if (!length(levels(outcome)) == 2L) 
    stop(sprintf("outcome %s must be binary", outcome.field))
 # cat(sprintf("\"%s\" will be the favorable outcome \n", levels(outcome)[2]))
  d[, outcome.field.which] <- outcome
  treatment <- d[, treatment.field.which]
  if (!is.numeric(treatment) & !is.integer(treatment)) 
    treatment <- as.numeric(treatment)
  if (!(length(unique(treatment)) == 2L & all(c(0, 1) %in% 
                                              unique(treatment)))) 
    stop(sprintf("%s, response must be numeric:\n 0 = no treatment \n 1 = treatment \n", 
                 treatment.field))
  d[, treatment.field.which] <- treatment
  predictors <- colnames(dataset)[-c(outcome.field.which, treatment.field.which)]
  predictors.next <- vector()
  iter = 1
  for (i in predictors) {
    iter <- length(predictors.next) + 1
    var <- d[, i]
    if (is.numeric(var) | is.integer(var)) {
      predictors.next[iter] <- i
    }
    if (is.character(var)) {
      var <- as.factor(var)
    }
    if (is.factor(var)) {
      if (length(levels(var)) > 2) {
        if (isTRUE(interactions)) {
          cat(sprintf("Creation of dummy variables for %s \n", 
                      i))
          for (l in levels(var)) {
            n <- paste(i, l, sep = "_")
            d[, n] <- ifelse(var == l, 1, 0)
            predictors.next[iter] <- n
            cat(sprintf("Dummy variable %s created \n", 
                        n))
            iter <- iter + 1
          }
        }
        else {
          if (length(levels(var)) > 32) {
            stop(cat(sprintf("%s has too many levels (superior to 32)", 
                             i)))
          }
          else {
            warning(sprintf("%s has more than 2 levels. Virtual Twins won't be able to run with interactions.", 
                            i))
            predictors.next[iter] <- i
          }
        }
      }
      else if (length(levels(var)) == 2) {
        cat(sprintf("%s is two-level factor. It has to be transformed into numeric value : \n", 
                    i))
        cat(sprintf("%s becomes 0 \n", levels(var)[1]))
        cat(sprintf("%s becomes 1 \n", levels(var)[2]))
        d[, i] <- ifelse(var == levels(var)[1], 0, 1)
        predictors.next[iter] <- i
      }
      else {
        cat(sprintf("%s is deleted because only one level", 
                    i))
      }
    }
  }
  colnames.order <- c(outcome.field, treatment.field, predictors.next)
  d <- d[, colnames.order]
  return(invisible(d))
}

  vt.data2<-function (dataset, outcome.field, treatment.field, interactions = TRUE, 
          ...) 
{
  data <- formatRCTDataset2(dataset, outcome.field, treatment.field, 
                           interactions = TRUE)
  VT.object(data = data, ...)
  }
  

  
  
  # Lebesque-Stieljes integral
  LS.int.FG<-function(fatx,x){
    temp<-fatx[1:(length(fatx)-1)]
    d.x<-diff(x)
    cum.f<-x[1]+cumsum(d.x*temp)
    return(c(x[1],cum.f))
  }
  
  
  
  FG.transform<-function(time,censor,S.cens,g.gamma=0.0,tail=0.0){
    n.censored<-sum(censor==1)
    if(n.censored<=3){
      warning("Number of censored observations less than 3, no transformation applied (return original times)")
      g.time<-time
    }
    else{
      # Calculate Integral
      # First modify S.cens
      S.mod<-ifelse(S.cens>g.gamma,S.cens,g.gamma)
      integrand.1<-1/S.mod
      
      int.1<-LS.int.FG(integrand.1,time)
      
      num.theta<-int.1-time
      den.theta<-(time/S.mod)-int.1
      
      temp<-num.theta[censor==1]/den.theta[censor==1]
      theta<-min(temp)
      
      temp2<-ifelse(censor==1,(1+theta)*int.1,(1+theta)*int.1-(theta*time/S.mod))
      g.time<-ifelse(S.cens>=tail,temp2,time)
    }
    return(g.time)
  }
  
  get.FG<-function(df,tte.name,event.name){
    df$id<-c(1:dim(df)[1])
    df.new<-df
    df1<-df[order(df[,tte.name]),]
    id1<-df1$id
    tt<-df1[,tte.name]
    cc<-1-df1[,event.name]
    # K-M of censoring distribution (Events are now censorings)
    #temp<-NA.CHR(X=tt,D=cc)
    temp<-NA.CHR.Weighted(time=tt,Delta=cc)
    g.time<-FG.transform(time=tt,censor=cc,S.cens=temp$S.NA,g.gamma=0.1,tail=0.1)
    id.match<-match(df$id,id1)
    df.new$Y.FG<-g.time[id.match]
    df.new$event.YFG<-rep(1,length(g.time))
    return(df.new)
  }

  
# vt.subg.harm.survival: Virtual twins approach for identifying
# 'harm' subgroup based on binary endpoint of survival
# cf.names = confounder (baseline) names
# y.name = outcome name
# treat.name = treatment (0/1) name
  
vt.subg.harm.survival<-function(data,cf.names,tte.name,event.name,treat.name,seed=8316951,data.pred=NULL,
                                details=TRUE,n.min=80,treat.threshold=0.0,maxdepth=3,tau=4,id.harm=NULL,ntree=500){

set.seed(seed)

treat<-data[,treat.name]
y<-data[,tte.name]
event<-data[,event.name]

y.tau<-as.numeric(ifelse(y<=tau,1,0))
y.tau[which(event==0 & y<=tau)]<-NA
data$y.tau<-y.tau

if(details){
  ndelete<-sum(is.na(y.tau))
  if(ndelete>0){
  cat("# of excluded censored patients due to binary endpoint=",c(ndelete),"\n")
  }
}

# Convert factors to numeric
xmat<-as.matrix(data[,cf.names])
xmat2<-NULL
for(kk in 1:ncol(xmat)){
  xmat2<-cbind(xmat2,c(as.numeric(xmat[,kk])))
}
colnames(xmat2)<-colnames(xmat)

#xmat<-as.matrix(data[,cf.names])
#xmat2<-apply(xmat,2,as.numeric)

vt.gsub<-NULL
diff.max<-NULL
sg.diff<-NULL
if(is.null(id.harm)){

# Exclude any missing data
df<-na.omit(data.frame(y.tau,treat,xmat2))

# vt.data2 is same as vt.data except automatic printing is removed
vt.o<-vt.data2(df, "y.tau", "treat", TRUE)

vt.f.rf <- vt.forest("one", vt.data = vt.o, interactions = T, ntree = ntree)

# grow RF for T = 1
model.rf.trt1 <- randomForest(x = vt.o$getX(trt = 1),y = vt.o$getY(trt = 1))
# grow RF for T = 0
model.rf.trt0 <- randomForest(x = vt.o$getX(trt = 0),y = vt.o$getY(trt = 0))
# initialize VT.forest.double()
vt.doublef.rf <- vt.forest("double",vt.data = vt.o,model_trt1 = model.rf.trt1, model_trt0 = model.rf.trt0)

twin1_random <- runif(length(y))
twin2_random <- runif(length(y))

model.difft <- VT.difft(vt.o, twin1 = twin1_random, twin2 = twin2_random, "absolute")
model.difft$computeDifft()

tr.class <- vt.tree("class",vt.difft = vt.f.rf,
                    sens = ">",
                    threshold = quantile(vt.f.rf$difft, seq(.5, .8, .1)),
                    maxdepth = maxdepth,cp = 0,maxcompete = 2)


vt.gsub <- vt.subgroups(tr.class)

# Minimum group size
group.size<-as.numeric(vt.gsub[,2])
treat.rate<-as.numeric(vt.gsub[,3])

vt.gsub<-vt.gsub[which(group.size>=n.min & treat.rate>=treat.threshold),]
# Note: Treatement (typo) is per the package
sg.diff<-as.numeric(c(vt.gsub$"Treatement event rate"))-as.numeric(c(vt.gsub$"Control event rate"))
# Order by sg.diff
vt.gsub<-vt.gsub[order(sg.diff),]

sg.diff<-as.numeric(c(vt.gsub$"Treatement event rate"))-as.numeric(c(vt.gsub$"Control event rate"))

if(details){ 
print(vt.gsub)
print(sg.diff)
}

loc.harm<-which.max(sg.diff)

id.harm<-vt.gsub$Subgroup[loc.harm]
diff.max<-sg.diff[loc.harm]
}

# Original analysis dataset is output with harm flag
data$treat.recommend<-NA
data$ID<-c(1:nrow(data))
# xmat2 to identify subgroups
# Add harm flag to original dataset
temp<-data.frame(xmat2)
temp$ID<-c(1:nrow(temp))
df.harm<-subset(temp,eval(parse(text=id.harm)))

pt.harm<-data$ID %in% df.harm$ID
data[pt.harm,"treat.recommend"]<-0

id.noharm<-setdiff(temp$ID,df.harm$ID)
pt.noharm<-data$ID %in% temp[id.noharm,]$ID
data[pt.noharm,"treat.recommend"]<-1

# For prediction dataset (e.g., "cross-validation testing dataset")
if(!is.null(data.pred)){
  
  treat<-data.pred[,treat.name]
  y<-data.pred[,tte.name]
  event<-data.pred[,event.name]
  
  y.tau<-as.numeric(ifelse(y<=tau,1,0))
  y.tau[which(event==0 & y<=tau)]<-NA
  data.pred$y.tau<-y.tau

  # Convert factors to numeric
  xmat<-as.matrix(data.pred[,cf.names])
  xmat3<-NULL
  for(kk in 1:ncol(xmat)){
    xmat3<-cbind(xmat3,c(as.numeric(xmat[,kk])))
  }
  colnames(xmat3)<-colnames(xmat)
  
  #xmat<-as.matrix(data.pred[,cf.names])
  #xmat3<-apply(xmat,2,as.numeric)
  
  data.pred$treat.recommend<-NA
  data.pred$ID<-c(1:nrow(data.pred))
  # xmat3 to identify subgroups
  # Add harm flag to original dataset
  temp<-data.frame(xmat3)
  temp$ID<-c(1:nrow(temp))
  df.harm.pred<-subset(temp,eval(parse(text=id.harm)))
  
  pt.harm.pred<-data.pred$ID %in% df.harm.pred$ID
  data.pred[pt.harm.pred,"treat.recommend"]<-0
  
  id.noharm.pred<-setdiff(temp$ID,df.harm.pred$ID)
  pt.noharm.pred<-data.pred$ID %in% temp[id.noharm.pred,]$ID
  data.pred[pt.noharm.pred,"treat.recommend"]<-1
}  
  
return(list(data=data,data.pred=data.pred,vt.gsub=vt.gsub,sg.harm.factors=id.harm,diff.max=diff.max,sg.diff=sg.diff))
}



vt.estimates.out<-function(df,vt.est=NULL,dgm=NULL,cox.formula.sim=NULL,cox.formula.adj.sim=NULL,vt.threshold=0.1,analysis="VT"){
  p.cens<-mean(1-df[,event.name])  
  
  Y<-df[,outcome.name]
  W<-df[,treat.name]
  D<-df[,event.name]
  taumax<-min(c(max(Y[W==1 & D==1]),max(Y[W==0 & D==1])))
  
  
  if(!is.null(vt.est)){
  vt.diff<-vt.est$diff.max
  sg.harm.vt<-NULL
  if(vt.diff>=vt.threshold) sg.harm.vt<-vt.est$sg.harm.factors
  }
  
  if(is.null(vt.est)){
    vt.diff<-NA
    sg.harm.vt<-NULL
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


  if(dgm$model=="null" & !is.null(dgm$vt.harm.true)) stop("For dgm model null vt.harm.true should be null")
  
  # dgm$fs.harm.true=NULL --> H does not exist
  
  # If H exists and found something
  if(!is.null(dgm$vt.harm.true) & !is.null(sg.harm.vt)){
    any.H<-1.0
    dfout<-vt.est$data
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
    #hr.H.true<-c(summary(coxph(cox.formula.sim,data=subset(df,flag.harm==1),robust=TRUE))$conf.int[1])
    #hr.Hc.true<-c(summary(coxph(cox.formula.sim,data=subset(df,flag.harm==0),robust=TRUE))$conf.int[1])
    
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
    #hr.H.hat<-c(summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==0),robust=TRUE))$conf.int[1])
    #hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==1),robust=TRUE))$conf.int[1])
    
    
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
    found.1<-sum(grepl(dgm$vt.harm.true[1],sg.harm.vt))
    found.2<-sum(grepl(dgm$vt.harm.true[2],sg.harm.vt))
    found.both<-as.numeric(c(found.1+found.2)==2)
    found.al3<-as.numeric(length(sg.harm.vt)>=3)
    
  }
  
  # If H exists and did NOT find something
  if(!is.null(dgm$vt.harm.true) & is.null(sg.harm.vt)){
    any.H<-0
    size.H<-0
    size.Hc<-nrow(df)
    
  
    # hr estimate with true H
    #hr.H.true<-c(summary(coxph(cox.formula.sim,data=subset(df,flag.harm==1),robust=TRUE))$conf.int[1])
    #hr.Hc.true<-c(summary(coxph(cox.formula.sim,data=subset(df,flag.harm==0),robust=TRUE))$conf.int[1])
    
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
    #hr.H.hat<-NA
    # Define HC as ITT 
    #hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=df,robust=TRUE))$conf.int[1])
    
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
  if(is.null(dgm$vt.harm.true) & !is.null(sg.harm.vt)){
    any.H<-1.0
    dfout<-vt.est$data
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
    #hr.H.true<-NA
    # ITT
    #hr.Hc.true<-c(summary(coxph(cox.formula.sim,data=df,robust=TRUE))$conf.int[1])
    
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
    #hr.H.hat<-c(summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==0),robust=TRUE))$conf.int[1])
    #hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==1),robust=TRUE))$conf.int[1])
    
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
    found.1<-sum(grepl(c("v1>=0.5"),sg.harm.vt))
    found.2<-sum(grepl(c("v4>=0.5"),sg.harm.vt))
    found.both<-as.numeric(c(found.1+found.2)==2)
    found.al3<-as.numeric(length(sg.harm.vt)>=3)
    
  }
  
  # If H does NOT exist and did NOT find something
  # Define hat(Hc)=ITT
  if(is.null(dgm$vt.harm.true) & is.null(sg.harm.vt)){
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
    #hr.H.true<-NA
    # ITT
    #hr.Hc.true<-c(summary(coxph(cox.formula.sim,data=df,robust=TRUE))$conf.int[1])
    
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
    #hr.H.hat<-NA
    #hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=df,robust=TRUE))$conf.int[1])
    
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
    found.1<-sum(grepl(c("v1>=0.5"),sg.harm.vt))
    found.2<-sum(grepl(c("v4>=0.5"),sg.harm.vt))
    found.both<-as.numeric(c(found.1+found.2)==2)
    found.al3<-as.numeric(length(sg.harm.vt)>=3)
    
  }
  
  df.res<-data.frame(any.H,size.H,size.Hc,ppv,npv,specificity,sensitivity,found.1,found.2,found.both,found.al3,
                     hr.H.true,hr.Hc.true,hr.H.hat,hr.Hc.hat,b1.H,b2.H,b1.Hc,b2.Hc,p.cens,analysis,taumax,hr.itt,l.itt,u.itt,
                     hr.adj.itt,l.adj.itt,u.adj.itt,
                     l.H.true,u.H.true,
                     l.Hc.true,u.Hc.true,
                     l.H.hat,u.H.hat,
                     l.Hc.hat,u.Hc.hat)
  return(df.res)
}


# To-do

# For FG approach we use tte.name= OS, event.name= OS event
# However for VT analysis, the outcome is tte.name=Y.FG (event.name=event.YFG)
vt.subg.harm.survival.boot<-function(data,cf.names,tte.name,event.name,treat.name,tte.vt.name=tte.name,event.vt.name=event.name,details=TRUE,print.level=0,
                                     n.min=80,treat.threshold=0.0,maxdepth=3,tau=4,boots=20,cox.formula,vt.est,ntree=500,trim=0.0){
  t.start<-proc.time()[1]
  
 # harm factors estimate
 sg.harm.est<-vt.est$sg.harm.factors
  
 # Estimated hazard ratio for estimated subgroup= hat(H)
 # hat.theta(hat(H))
 df.0<-subset(vt.est$data,treat.recommend==0) # harm --> "recommend Control"
 df.1<-subset(vt.est$data,treat.recommend==1) # non-harm --> "recommend Treatment"
  
 coxfit.0<-summary(coxph(cox.formula,data=df.0,robust=TRUE))$conf.int
 coxfit.1<-summary(coxph(cox.formula,data=df.1,robust=TRUE))$conf.int
  
  # VT estimate (observed) for "recommend Control"
  hr.vt.0<-c(coxfit.0[1])
  # VT estimate (observed) for "recommend Treatment"
  hr.vt.1<-c(coxfit.1[1])
    # bootstrap data
  
  # Data for VT identification
  # Convert baseline factors to numeric
  xmat<-data.frame(data[,covs.model])
  xmat2<-xmat
  for(kk in 1:ncol(xmat2)){
    xmat2[,kk]<-as.numeric(xmat[,kk])
  }
  xmat2<-as.matrix(xmat2)
  
if(tte.vt.name!=tte.name)  df.vt<-data.frame(data[,c(tte.name,event.name,treat.name,tte.vt.name,event.vt.name)],xmat2)
if(tte.vt.name==tte.name)  df.vt<-data.frame(data[,c(tte.name,event.name,treat.name)],xmat2)
  
  # Bootstrap df.vt
  id0<-c(1:nrow(df.vt))
  df.vt$id<-id0
  
  bc.0.boots<-rep(NA,boots)
  bc.1.boots<-rep(NA,boots)
  
  vi.covs<-data.frame(confounders=covs.model,var.importance.boots=0)
  
  for(bb in 1:boots){
    
    set.seed(8316951+bb*1000)
    
    id.boot<-sample(id0,size=nrow(df.vt),replace=TRUE)
    
    df.vt.boot<-df.vt[id.boot,]
    
    df.vt.boot$ID<-c(1:nrow(df.vt.boot))
    
    # hatstar(theta)(hat(H)) ---> bootstrap hr evluated at hat(H)
    # Evaluate bootstrapped data at estimated H
    # id.harm=sg.harm.est (Setting H at estimate = hat(H))
    # data.pred = df.vt.boot (Evaluating boot data @ hat(H))
    
    temp<-try(vt.subg.harm.survival(data=data,cf.names=cf.names,tte.name=tte.vt.name,event.name=event.vt.name,treat.name=treat.name,
                                seed=8316951,data.pred=df.vt.boot,id.harm=sg.harm.est,details=details,
                                n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau,ntree=ntree),TRUE)
    
    if(!inherits(temp,"try-error")){ # Boot loop 1
      
    # Evaluate at data.pred
    df.0.boot<-subset(temp$data.pred,treat.recommend==0) # harm --> "recommend Control"
    df.1.boot<-subset(temp$data.pred,treat.recommend==1) # non-harm --> "recommend Treatment"
    
    coxfit.0.boot<-summary(coxph(cox.formula,data=df.0.boot,robust=TRUE))$conf.int
    coxfit.1.boot<-summary(coxph(cox.formula,data=df.1.boot,robust=TRUE))$conf.int
    
    # VT estimate (observed) for "recommend Control"
    hr.vt.0.boot<-c(coxfit.0.boot[1])
    # VT estimate (observed) for "recommend Treatment"
    hr.vt.1.boot<-c(coxfit.1.boot[1])
    
    # Bias correction term1
    bc.0.term1<-hr.vt.0.boot-hr.vt.0
    bc.1.term1<-hr.vt.1.boot-hr.vt.1
    
    # Estimate H based on boostrap data
    # Want to evaluate based on observed data
    # Set data.pred = df.vt
    temp<-try(vt.subg.harm.survival(data=df.vt.boot,cf.names=cf.names,tte.name=tte.vt.name,event.name=event.vt.name,treat.name=treat.name,
                                seed=8316951,data.pred=df.vt,details=details,
                                n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau,ntree=ntree),TRUE)
    
    if(!inherits(temp,"try-error")){ # Boot loop 2
      
    # hatstar(H)
    sg.harm.boot<-temp$sg.harm.factors
    
    for(jj in 1:length(cf.names)){ 
      vi.covs[jj,2]<-vi.covs[jj,2]+as.numeric(grepl(cf.names[jj],sg.harm.boot))
    }
    
    # Evaluate at data.pred=df.vt
    df.0.boot.hstar<-subset(temp$data.pred,treat.recommend==0) # harm --> "recommend Control"
    df.1.boot.hstar<-subset(temp$data.pred,treat.recommend==1) # non-harm --> "recommend Treatment"
    
    coxfit.0.boot.hstar<-summary(coxph(cox.formula,data=df.0.boot.hstar,robust=TRUE))$conf.int
    coxfit.1.boot.hstar<-summary(coxph(cox.formula,data=df.1.boot.hstar,robust=TRUE))$conf.int
    
    # VT estimate (observed) for "recommend Control"
    hr.vt.0.boot.hstar<-c(coxfit.0.boot.hstar[1])
    # VT estimate (observed) for "recommend Treatment"
    hr.vt.1.boot.hstar<-c(coxfit.1.boot.hstar[1])
    
    
    # Bias correction term2
    bc.0.term2<-hr.vt.0.boot.hstar-hr.vt.0
    bc.1.term2<-hr.vt.1.boot.hstar-hr.vt.1
    
    bc.0.boots[bb]<-bc.0.term1+bc.0.term2
    bc.1.boots[bb]<-bc.1.term1+bc.1.term2
  }
  
  } # End boot loop 1
  } # End boot loop 2
  
  bc.0.boots<-na.omit(bc.0.boots)
  bc.1.boots<-na.omit(bc.1.boots)
  
  if(trim>0){
  bc.0.boots<-x.truncate(bc.0.boots,truncate=trim)
  bc.1.boots<-x.truncate(bc.1.boots,truncate=trim)
  }  
  
  vi.covs<-na.omit(vi.covs)
  n.boots<-length(bc.0.boots)
  
  hr.vt.0.bc<-hr.vt.0-mean(bc.0.boots)
  hr.vt.1.bc<-hr.vt.1-mean(bc.1.boots)
  
  if(print.level>0){
    t.end<-proc.time()[1]
    t.min<-(t.end-t.start)/60
    
    cat("# of bootstrap estimation samples, timing (minutes)",c(n.boots,t.min),"\n")
    cat("Observed estimate: Rec. Control, Rec. Treat",c(hr.vt.0,hr.vt.1),"\n")
    cat("Bias corrected estimate: Rec. Control, Rec. Treat",c(hr.vt.0.bc,hr.vt.1.bc),"\n")
  }
  
  vi.covs<-vi.covs[order(vi.covs[,2],decreasing=TRUE),]
  vi.covs[,2]<-vi.covs[,2]/n.boots
  
  
  return(list(hr.vt.0.bc=hr.vt.0.bc,bc.0.boots=bc.0.boots,vi.covs=vi.covs,
  hr.vt.1.bc=hr.vt.1.bc,bc.1.boots=bc.1.boots,hr.vt.0.boots=c(rep(hr.vt.0,n.boots)-bc.0.boots),
  hr.vt.1.boots=c(rep(hr.vt.1,n.boots)-bc.1.boots)))
  
}










vt.subg.harm.binary.combine<-function(data,cf.names,y.name,treat.name,seed=8316951,data.pred=NULL,c.threshold=0.05,ntree=500,
                                details=TRUE,n.min=80,treat.threshold=0.0,maxdepth=3,id.harm=NULL){
  
  set.seed(seed)
  
  treat<-data[,treat.name]
  y<-data[,y.name]
  
  # Convert factors to numeric
  xmat<-data.frame(data[,covs.model])
  xmat2<-xmat
  for(kk in 1:ncol(xmat2)){
    xmat2[,kk]<-as.numeric(xmat[,kk])
  }
  xmat2<-as.matrix(xmat2)
  
  vt.gsub<-NULL
  diff.max<-NULL
  sg.diff<-NULL
  vt.gsub2<-NULL
  if(is.null(id.harm)){
    
    # Exclude any missing data
    df<-na.omit(data.frame(y,treat,xmat2))
    
    # vt.data2 is same as vt.data except automatic printing is removed
    vt.o<-vt.data2(df, "y", "treat", TRUE)
    
    vt.f.rf <- vt.forest("one", vt.data = vt.o, interactions = T, ntree = ntree)
    
    # grow RF for T = 1
    model.rf.trt1 <- randomForest(x = vt.o$getX(trt = 1),y = vt.o$getY(trt = 1))
    # grow RF for T = 0
    model.rf.trt0 <- randomForest(x = vt.o$getX(trt = 0),y = vt.o$getY(trt = 0))
    # initialize VT.forest.double()
    vt.doublef.rf <- vt.forest("double",vt.data = vt.o,model_trt1 = model.rf.trt1, model_trt0 = model.rf.trt0)
    
    twin1_random <- runif(length(y))
    twin2_random <- runif(length(y))
    
    model.difft <- VT.difft(vt.o, twin1 = twin1_random, twin2 = twin2_random, "absolute")
    model.difft$computeDifft()
    
    tr.class <- vt.tree("class",vt.difft = vt.f.rf,
                        sens = ">",
                        threshold = quantile(vt.f.rf$difft, seq(.5, .8, .1)),
                        maxdepth = maxdepth,cp = 0,maxcompete = 2)
    
    
    vt.gsub <- vt.subgroups(tr.class)
    
    #print(vt.gsub)
    
    # Minimum group size
    group.size<-as.numeric(vt.gsub[,2])
    treat.rate<-as.numeric(vt.gsub[,3])
    
    # Find at least one (reset n.min to max(group.size))
    if(max(group.size)<n.min) n.min<-max(group.size)
    
    vt.gsub<-vt.gsub[which(group.size>=n.min & treat.rate>=treat.threshold),]
    # Note: Treatement (typo) is per the package
    
    if(!is.null(vt.gsub)){
    sg.diff<-as.numeric(c(vt.gsub$"Treatement event rate"))-as.numeric(c(vt.gsub$"Control event rate"))
    # Order by sg.diff
    vt.gsub<-vt.gsub[order(sg.diff),]
    
    sg.diff<-as.numeric(c(vt.gsub$"Treatement event rate"))-as.numeric(c(vt.gsub$"Control event rate"))
    
    if(details){ 
      print(vt.gsub)
      print(sg.diff)
    }
   
    loc.harm<-which.max(sg.diff)
    id.harm<-vt.gsub$Subgroup[loc.harm]
    diff.max<-sg.diff[loc.harm]
    
    # If no subgroups satisfy c.threshold --> chose only max  
    if(diff.max<c.threshold) c.threshold<-diff.max
    
    # Among all subgroups with sg.diff>=c.threshold 
    vt.gsub2<-vt.gsub[sg.diff>=c.threshold,]
    
  if(details)  print(vt.gsub2)
    
    }
  }
  
  if(!is.null(vt.gsub2)){
  # Original analysis dataset is output with harm flag
  data$treat.recommend<-NA
  data$ID<-c(1:nrow(data))
  # xmat2 to identify subgroups
  # Add harm flag to original dataset
  temp<-data.frame(xmat2)
  temp$ID<-c(1:nrow(temp))
  
  df.harm<-NULL
  for(ii in 1:nrow(vt.gsub2)){
  id.harmi<-vt.gsub2$Subgroup[ii]
  #print(id.harmi)
  df.harm<-rbind(df.harm,subset(temp,eval(parse(text=id.harmi))))
                 #print(nrow(df.harm))
  }               
  
  pt.harm<-data$ID %in% df.harm$ID
  data[pt.harm,"treat.recommend"]<-0
  
  id.noharm<-setdiff(temp$ID,df.harm$ID)
  pt.noharm<-data$ID %in% temp[id.noharm,]$ID
  data[pt.noharm,"treat.recommend"]<-1
  
  # For prediction dataset (e.g., "testing dataset")
  if(!is.null(data.pred)){
    
    treat<-data.pred[,treat.name]
    y<-data.pred[,y.name]

    # create flag for prediction dataset
    # Convert factors to numeric
    xmat<-data.frame(data.pred[,covs.model])
    xmat3<-xmat
    for(kk in 1:ncol(xmat3)){
      xmat3[,kk]<-as.numeric(xmat[,kk])
    }
    xmat3<-as.matrix(xmat3)
    
    data.pred$treat.recommend<-NA
    data.pred$ID<-c(1:nrow(data.pred))
    # xmat3 to identify subgroups
    # Add harm flag to original dataset
    temp<-data.frame(xmat3)
    temp$ID<-c(1:nrow(temp))
    
    df.harm.pred<-NULL
    for(ii in 1:nrow(vt.gsub2)){
      id.harmi<-vt.gsub2$Subgroup[ii]
     # print(id.harmi)
            df.harm.pred<-rbind(df.harm.pred,subset(temp,eval(parse(text=id.harmi))))
          }               
    pt.harm.pred<-data.pred$ID %in% df.harm.pred$ID
    data.pred[pt.harm.pred,"treat.recommend"]<-0
    
    id.noharm.pred<-setdiff(temp$ID,df.harm.pred$ID)
    pt.noharm.pred<-data.pred$ID %in% temp[id.noharm.pred,]$ID
    data.pred[pt.noharm.pred,"treat.recommend"]<-1
  }  
  }
  return(list(data=data,data.pred=data.pred,vt.gsub=vt.gsub,sg.harm.factors=id.harm,diff.max=diff.max,sg.diff=sg.diff))
}


