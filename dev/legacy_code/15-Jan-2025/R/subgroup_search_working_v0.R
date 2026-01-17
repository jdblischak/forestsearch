# minp = Minimum prevalence rate 
# Only combinations where all factors have prevalance at least minp are evaluated
# Minimum difference in subgroup sample size
# Max subgroup size is n: After first factor x1 is added, the sample size for 
# the subgroup is then n(x1), say.   If combining with next factor x2 
# does not reduce the sample size by rmin, then consider this combination "redundant"
# Adding d.min to require min(Events) per treatment arm (d0.min for control, d1.min for treatment)
# covs.in[isin,]<-unlist(lapply(covs.in[isin,],function(x){ifelse(x==1,0,1)})) # slower

# Two-factor combinations
# We just need a placeholder for possible combinations 
# where at most 2 factor combinations are evaluated
#
#x1<-seq(0,1)  # Possible binomial realizations from sample of size n            
#L <- 10 (5 factors, 2 levels each)
#xx<-expand.grid(x1,x1,x1,x1,x1,x1,x1,x1,x1,x1)
#xx<-expand.grid(x1,x1,x1,x1,x1,x1,x1,x1)
#tot <- apply(xx,1,sum)
#sum(tot >0 & tot <= 2)
# L = 8 -- > 36
#(2**L)-1
# For L=10 there are 55
# L=10 -> 55
# L=12 --> 78
# L=13 --> 91
# L=20 --> 210 (10 factors, 2 levels each)
# L=30 (15 factors 2 levels each)
# Not feasible in this manner
# Set placeholder to 2,000

# Removing ID storage since not needed
subgroup.search<-function(Y,Event,Treat,ID=NULL,Z,n.min=30,d0.min=15,d1.min=15,
hr.threshold=1.0,max.minutes=30,minp=0.05,rmin=5,details=FALSE,maxk=2,min_screen_target=NULL){

non.redundant<-0.0
ngroups.fit<-0.0
# Number of SG's meeting overall criteria with estimable hr
ngroups.meet<-0.0
# Minimal target for ngroups.meet = min_screen_target
  
  # output max of sg estimates
  # initialize
  max_sg_est<--Inf
  
  id.group<-0.0 # Identify subgroup (Will output subgroups directly)
  subs.out<-NULL
  
  temp<-cbind(Y,Event,Treat,Z)
  temp<-na.exclude(temp)
  yy<-temp[,1]
  dd<-temp[,2]
  tt<-temp[,3]
  zz<-temp[,-c(1,2,3)]
  
  n<-length(yy)
  L<-ncol(zz)
  covs.in<-as.matrix(rep(0,L)) # Initially no members
  # Store results appending to HR.model.k
  HR.model.k<-NULL

  limit_all<-(2**L)-1 # All-possible configurations including redundancies 
  
  # Maximum number of combinations among factors <= maxk
  if(maxk==1) max_count <- L 
  if(maxk==2){
  # choose(L,2)+choose(L,1)  
  max_count <- (L*(L-1)/2) + L  
  }
  if(maxk==3){
  choose(L,3)+choose(L,2)+choose(L,1)
  max_count <- ((L*(L-2)*(L-1))/6)+(L*(L-1)/2) + L  
  }
  
 if(details){
    cat("Number of unique levels (L) and all-possible subgroups=",c(L,limit_all),"\n")
    if(limit_all>10^6) cat("Number of all-possible subgroups (in millions)=",c(limit_all/(10^6)),"\n")
   cat("Number of possible configurations <= maxk: maxk, #=",c(maxk,max_count),"\n")
   }
  # Track criteria not met for subgroup evaluation
  crit0.met <- 0.0
  crit.fail0<-0
  # crit.fail0 is first criteria to be met
  # The following criteria are conditional on meeting crit.fail0
  crit.failure<-0
  crit.faild0<-0
  crit.faild1<-0
  crit.failn<-0
  t.start<-proc.time()[3]
  t.sofar<-0
  k_count <- 0.0
  
  index_2factor <- t(combn(L, 2))
  counts_2factor <- nrow(index_2factor)
  index_1factor <- t(combn(L, 1))
  counts_1factor <- nrow(index_1factor)
  
  if(!is.null(min_screen_target)) min_target <- min_screen_target
  if(is.null(min_screen_target)) min_target <- max_count
  
  tot_counts <- counts_2factor+counts_1factor
  if(tot_counts != max_count) stop("Error with maximum combinations kmax=2")
  
  # Loop through single-factor subgroups first
  
  for(kk in 1:tot_counts){
  covs.in<-as.matrix(rep(0,L)) # Initially no members
  if(kk <= counts_1factor){
  which1 <- index_1factor[kk]  
  covs.in[which1] <- 1.0
  }
  if(kk > counts_1factor){
  # reset kindex
  kk_new <- (kk-counts_1factor) 
  which1 <- index_2factor[kk_new,1]
  which2 <- index_2factor[kk_new,2]
  covs.in[c(which1,which2)] <- c(1.0,1.0)
  }
    ## We now have subgroup indices
    ## So extract subgroup membership
    k.in<-sum(covs.in)  
    # check timing
    t.now<-proc.time()[3]
    t.sofar<-c((t.now-t.start)/60)
    # Stop search algorithm if time exceeded: Output accumulated results
    if(t.sofar>max.minutes){ 
    out<-NULL
   if(details){
    cat("Time exceeding max.minutes, stopping search algorithm and continuing with SG's found so far [if any]","\n")
    }   
      if(!is.null(HR.model.k)){
        hr.out<-data.table(HR.model.k)
        names(hr.out)<-c("grp","K","n","E","d1","m1","m0","HR","L(HR)","U(HR)",colnames(Z))
        rownames(hr.out)<-NULL
        hr.out<-setorder(hr.out,-HR,K)
        out.found<-list(ngroups.fit=ngroups.fit,hr.subgroups=hr.out)
        out<-list(out.found=out.found,time_search=t.sofar,
                  max_sg_est=c(hr.out[1,"HR"]),L=L)
      }
      if(details){ 
        cat("# of subgroups based on 2-factor combinations",c(crit0.met),"\n")
        cat("% of all-possible combinations (<= maxk)",c(100*round(crit0.met/max_count,3)),"\n")
        cat("# of subgroups=",c(non.redundant),"\n")
        cat("# of subgroups fitted=",c(ngroups.fit),"\n")
        cat("Minutes=",c(t.sofar),"\n")
      }    
      return(out)
      stop("Maximum timing exceeded")
    }
    crit.fail0 <- crit.fail0+(k.in>maxk)
    
    if(k.in <= maxk){
    crit0.met <- crit0.met+1
      
    if(details){
      
      if(crit0.met==ceiling(max_count/20)){
        cat("Approximately 5% of max_count met: minutes",c(t.sofar),"\n")
      }
      if(crit0.met==ceiling(max_count/10)){
        cat("Approximately 10% of max_count met: minutes",c(t.sofar),"\n")
      }
      if(crit0.met==ceiling(max_count/5)){
        cat("Approximately 20% of max_count met: minutes",c(t.sofar),"\n")
      }
    if(crit0.met==ceiling(max_count/3)){
      cat("Approximately 33% of max_count met: minutes",c(t.sofar),"\n")
    }
    if(crit0.met==ceiling(max_count/2)){
      cat("Approximately 50% of max_count met: minutes",c(t.sofar),"\n")
    }
      if(crit0.met==ceiling(3*max_count/4)){
        cat("Approximately 75% of max_count met: minutes",c(t.sofar),"\n")
      }  
      if(crit0.met==ceiling(0.9*max_count)){
        cat("Approximately 90% of max_count met: minutes",c(t.sofar),"\n")
      }  
    }  
    
    x<-matrix(NA,n,k.in)
    index<-1
      for(jj in 1:L){
        if(covs.in[jj]==1){
          x[,index]<-zz[,jj]
          index<-index+1
        } 
      }
      # If any two factors are mutually exclusive then set id.x=NULL (-> next iteration)  
      xpx<-t(x)%*%x
      if(!all(xpx>0)) id.x<-NULL
       if(all(xpx>0)){
        id.x<-rep(1,n)
        # Elements will be 1 if in subgroup defined by x  
        flag.redundant<-FALSE
        nx.last <- n # To flag redundancy
        for(m in 1:ncol(x)){
          if(!flag.redundant){
           id.x <- c(id.x)*c(x[,m])
            nx.this <- sum(id.x)
           if(nx.last-nx.this <= rmin) flag.redundant<-TRUE
            nx.last <- nx.this
         }
        } # flag redundant
      }  
      
      d1 <- sum(c(dd*tt)[which(id.x==1)]) # number of events for treatment
      d0 <- sum(c(dd*(1-tt))[which(id.x==1)]) # number of events for control
      nx<-sum(id.x)
      crit.failn<- crit.failn +c(nx<=n.min)
      crit.faild1<-crit.faild1+c(d1<d1.min)
      crit.faild0<-crit.faild0+c(d0<d0.min)
      crit3<-(all(apply(x,2,mean)>=minp))
      
      if(!flag.redundant & d1>=d1.min & d0>=d0.min & nx>n.min & crit3){
      data.x<-data.table(Y=yy,E=dd,Treat=tt,id.x=id.x)
        
       df.x<-subset(data.x,id.x==1)
       non.redundant<-non.redundant+1
    
        hr.cox<-try(summary(coxph(Surv(Y,E)~Treat,data=df.x,robust=FALSE))$conf.int,TRUE)
       
        if(!inherits(hr.cox,"try-error")){
         max_sg_est<-max(c(max_sg_est,hr.cox[1]))
         ngroups.fit<-ngroups.fit+1
          
          if(hr.cox[1]>hr.threshold){
            ngroups.meet<-ngroups.meet+1
            meds<-summary(survfit(Surv(Y,E)~Treat,data=df.x))$table[,"median"]
            m0<-meds[1]
            m1<-meds[2]
            id.group<-id.group+1
            # ANY additions or modifications needs to update consistency function
            hr.this<-c(id.group,sum(covs.in),nx,d0+d1,d1,m1,m0,hr.cox[c(1,3,4)],covs.in)
            HR.model.k<-rbind(HR.model.k,hr.this)
          } # Estimation  criterion loop
        }
      }
  } # maxk 
  } # d1, d0, nmin | maxk criteria
  t.end<-proc.time()[3]
  t.min<-(t.end-t.start)/60
  
  prop_max_count=c(100*round(crit0.met/max_count,3))
  
  if(details){
    cat("k_count",c(k_count),"\n")
    cat("# of subgroups based on 2-factor combinations",c(crit0.met),"\n")
    cat("% of all-possible combinations (<= maxk)",c(100*round(crit0.met/max_count,3)),"\n")
    cat("# of subgroups based on # variables > k.max and excluded (per million)",c(sum(crit.fail0)/10^6),"\n")
    cat("k.max=",c(maxk),"\n")
    cat("Events criteria for control,exp=",c(d0.min,d1.min),"\n")
    cat("# of subgroups with events less than criteria: control, experimental",c(sum(crit.faild0),sum(crit.faild1)),"\n")
    cat("# of subgroups with sample size less than criteria",c(sum(crit.failn)),"\n")
    cat("# of subgroups meeting all criteria =",c(non.redundant),"\n")
    cat("# of subgroups fitted (Cox model estimable) =",c(ngroups.fit),"\n")
    cat("*Subgroup Searching Minutes=*",c(t.min),"\n")
    cat("Number of subgroups meeting HR threshold",c(ngroups.meet),"\n")
    }
  
  out.found<-NULL
  if(!is.null(HR.model.k)){
    hr.out<-data.table(HR.model.k)
    names(hr.out)<-c("grp","K","n","E","d1","m1","m0","HR","L(HR)","U(HR)",colnames(Z))
    # ANY additions or modifications needs to update consistency function
    rownames(hr.out)<-NULL
    hr.out<-setorder(hr.out,-HR,K)
    
    out.found<-list(ngroups.fit=ngroups.fit,hr.subgroups=hr.out)
  }
  
  if(is.null(out.found) & details) cat("NO subgroup candidate found (FS)","\n")
  if(!is.null(out.found) & details) cat("Subgroup candidate(s) found (FS)","\n")
  
return(list(out.found=out.found,max_sg_est=max_sg_est,time_search=t.sofar,L=L,max_count=max_count,prop_max_count=prop_max_count))
}
