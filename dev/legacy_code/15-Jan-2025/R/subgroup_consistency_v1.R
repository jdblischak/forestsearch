
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

hr_splits<-function(nsplits,df,hrc){
foreach(bb=1:nsplits) %dorng% {  
# Random splitting
n.split<-round(nrow(df)/2,0)
    datax.id<-c(1:nrow(df))
    idx.split1<-sample(datax.id,size=n.split,replace=FALSE)
    df.x.split1<-df[idx.split1,]
    idx.split2<-setdiff(datax.id,idx.split1)
    df.x.split2<-df[idx.split2,]
    hr.split1<-try(suppressWarnings(summary(coxph(Surv(Y,Event)~Treat,data=df.x.split1,robust=FALSE))$conf.int),TRUE)
    hr.split2<-try(suppressWarnings(summary(coxph(Surv(Y,Event)~Treat,data=df.x.split2,robust=FALSE))$conf.int),TRUE)
    if(!inherits(hr.split1,"try-error") & !inherits(hr.split2,"try-error")){
    flag.consistency<-c(hr.split1[1,1]>hrc & hr.split2[1,1]> hrc)
    list(flag.consistency=flag.consistency)
    } # both splits estimable
  }
}


sg_consistency_out<-function(df,result_new,sg_focus,index.Z,names.Z,details=FALSE,plot.sg=FALSE,by.risk=4,confs_labels){

if(sg_focus=="hr") setorder(result_new,-Pcons,K)
# Choose sg with largest size 
# Emphasis by Nsg, consistency, and smallest K
if(sg_focus %in% c("Nsg","Nsg_only")) setorder(result_new,-N,-Pcons,K)
if(sg_focus %in% c("Msg","Msg_only")) setorder(result_new,N,-Pcons,K)
this.1<-NULL
data.found<-NULL
grp1<-NULL
  # show (up to) top 10
  if(details){
    cat("Number of subgroups meeting consistency criteria=","\n")
    print(result_new[1:min(10,nrow(result_new)),])
  }
  # Plot top 1
  # Subgroup corresponding to the top
  if(nrow(result_new)==1){
    top_result<-result_new  
  }
  if(nrow(result_new)>1){
    top_result<-result_new[1,]
  }
  grp1<-as.numeric(top_result$g)
  m1<-as.numeric(top_result$m)
  index1<-as.numeric(unlist(index.Z[m1,]))
  this.1<-names.Z[index1==1]
  
  this.1_label<-unlist(lapply(this.1,FS_labels,confs_labels=confs_labels))
  
  id.harm<-c(paste(this.1, collapse="==1 & "))
  id.harm<-paste(id.harm,"==1")
  df.sub<-subset(df,eval(parse(text=id.harm)))
  df.sub$treat.recommend<-0
  
  id.noharm<-c(paste(this.1, collapse="!=1 | "))
  id.noharm<-paste(id.noharm,"!=1")
  df.subC<-subset(df,eval(parse(text=id.noharm)))
  df.subC$treat.recommend<-1
  
  data.found<-rbind(df.sub,df.subC)
  
  # Original analysis dataset is output with harm flag
  
  if(details & plot.sg){
    plot.subgroup(sub1=df.sub,sub1C=df.subC,tte.name="Y",event.name="Event",treat.name="Treat",subid=paste(this.1_label,collapse=" & "),byrisk=by.risk)
    print(this.1_label)
    cat("% consistency criteria met=",c(top_result$Pcons),"\n")
  }
# For df.est --> only returd "id" and "treat.recommend"  
out<-list(result=result_new,sg.harm=this.1,sg.harm_label=this.1_label,df_flag=data.found[,c("id","treat.recommend")],sg.harm.id=grp1)
return(out)
}

# sg_focus = hr --> sort by HR's to prioritize largest HR's
# sg_focus= Nsg ---> prioritize SG size

# pstop_futile < 0 --> forces evaluation of all subgroups
# pstop_futile >=0 stops evaluation at first subgroup with consistency<=pstop_futile
# with idea that for sg_focus=HR the SG's are sorted by HR and thus once "futile"
# is likely to remain so for SG's with even lower HR's
# Note: is set to < 0 for sg_focus=Nsg or Msg in forestsearch()

# Remove sg_focus, will return all options
# Add Lsg = number of possible factors ("L" in paper)
# Add stop_Kgroups to stop when number of identified subgroups meeting consistency exceeds stop_Kgroups
# That is, Stopping for futility OR reaching stop_Kgroups (meeting Kgroups is combined with futility flag)
subgroup.consistency<-function(df,hr.subgroups,hr.threshold=1.0,hr.consistency=1.0,pconsistency.threshold=0.9,m1.threshold=Inf,n.splits=100,details=FALSE,
stop.threshold=1.1,by.risk=2,plot.sg=FALSE,maxk=7,pstop_futile=-0.1,Lsg,confs_labels,sg_focus="hr",stop_Kgroups=10){

names.Z<-c(names(hr.subgroups[,-c("grp","K","n","E","d1","m1","m0","HR","L(HR)","U(HR)")]))

if(length(names.Z)!=Lsg) stop("HR subgroup results not matching L, check subgroup search function")

if(m1.threshold==Inf) found.hrs<-hr.subgroups[hr.subgroups$HR>=hr.threshold,]  
  # Implement m1.threshold criteria for sg's where treatment medians are estimable
if(m1.threshold!=Inf){
hr.subgroups<-hr.subgroups[which(!is.na(hr.subgroups$m1)),]
found.hrs<-hr.subgroups[hr.subgroups$HR>=hr.threshold & hr.subgroups$m1<=m1.threshold,]
}

# sort by HRs unless Nsg_only or Msg_only
if(!(sg_focus %in% c("Nsg_only","Msg_only"))) found.hrs<-found.hrs[order(found.hrs$HR,decreasing=TRUE),]

if(sg_focus=="Nsg_only") found.hrs<-found.hrs[order(found.hrs$n,decreasing=TRUE),]
  
if(sg_focus=="Msg_only") found.hrs<-found.hrs[order(found.hrs$n,decreasing=FALSE),]

n.groups<-0 
  
# Extract sg factors
index.Z<-found.hrs[,-c("grp","K","n","E","d1","m1","m0","HR","L(HR)","U(HR)")]

if(dim(index.Z)[2]!=Lsg) stop("HR subgroup results not matching L, check subgroup search function")

result_new <- NULL
result <- NULL
this.1 <- NULL
# Store group id and consistency %
any.found <- 0
p.consistency <- 0
stop_futile <- FALSE

  if(details){
    t.start<-proc.time()[3]
    cat("SGs (1st 10) meeting screening thresholds sorted by sg_focus=",c(sg_focus),"\n")
    # Show up to first 6 factors
    nZ <- nrow(index.Z)
    n10 <- min(nZ,10)
    k10 <- min(nrow(found.hrs),10)
    temp <- found.hrs[,c("n","E","d1","HR","L(HR)")]
    temp <- as.data.frame(temp)
    temp <- cbind(temp,index.Z[,1:n10])
    print(round(temp[1:k10,],2))
  }  
  for(m in 1:nrow(found.hrs)){
  # Stopping at first sg exceeding stop.threshold
  # Or stopping after first sg below pstop_futile
  # Since HR's are sorted, chance of meeting 90% likely "futile"
    if(p.consistency < stop.threshold & !stop_futile){
      indexm<-as.numeric(unlist(index.Z[m,]))
      this.m<-names.Z[indexm==1]
      n.groups<-n.groups+1
      
      id.m<-c(paste(this.m, collapse="==1 & "))
      id.m<-paste(id.m,"==1")
      # data for subgroup "x"
      df.sub<-subset(df,eval(parse(text=id.m)))
      df.x<-data.table(df.sub)
      # Random splitting
      flag.consistency<-rep(0,n.splits) # Initialize as 0
      # If not estimable (both splits) then consider NOT consistent
      set.seed(8316951)
      for(bb in 1:n.splits){
        n.split<-round(nrow(df.x)/2,0)
        datax.id<-c(1:nrow(df.x))
        idx.split1<-sample(datax.id,size=n.split,replace=FALSE)
        df.x.split1<-df.x[idx.split1,]
        idx.split2<-setdiff(datax.id,idx.split1)
        df.x.split2<-df.x[idx.split2,]
       hr.split1<-try(suppressWarnings(summary(coxph(Surv(Y,Event)~Treat,data=df.x.split1,robust=FALSE))$conf.int),TRUE)
       hr.split2<-try(suppressWarnings(summary(coxph(Surv(Y,Event)~Treat,data=df.x.split2,robust=FALSE))$conf.int),TRUE)
        if(!inherits(hr.split1,"try-error") & !inherits(hr.split2,"try-error")){
        flag.consistency[bb]<-c(hr.split1[1,1]>hr.consistency & hr.split2[1,1]> hr.consistency)
          } # both splits estimable
       } # end splits 
      p.consistency<-mean(flag.consistency,na.rm=TRUE)
      
      # Trigger stopping for futility OR reaching stop_Kgroups
      if(p.consistency <= pstop_futile | (any.found==(stop_Kgroups-1))){
       stop_futile <- TRUE
      }
  
      if(details) cat("Consistency",p.consistency,"\n")
      
      if(p.consistency>pconsistency.threshold){
        
        any.found <- any.found+1

        k<-length(this.m)
        covsm<-rep("M",maxk) # store at most kmax names 
        mindex<-c(1:maxk)
        Mnames<-paste(covsm,mindex,sep=".")
        mfound<-matrix(rep("",maxk))
        
        this.m_label<-unlist(lapply(this.m,FS_labels,confs_labels=confs_labels))
      
        mfound[c(1:k)]<-this.m_label
        mfound<-c(mfound)
        
        #resultk<-c(p.consistency,found.hrs$n[m],found.hrs$grp[m],m,k,mfound)
        # Add events "E" to output 25-March-2025
        resultk<-c(p.consistency,found.hrs$HR[m],found.hrs$n[m],found.hrs$E[m],found.hrs$grp[m],m,k,mfound)
        
        resultk<-transpose(as.data.table(resultk))
        
        #names(resultk)<-c("Pcons","N","g","m","K",Mnames)
        
        names(resultk)<-c("Pcons","hr","N","E","g","m","K",Mnames)
        
        result<-rbind(result,resultk)
        
      if(details){
        cat("# of splits, Kmet=",c(n.splits,any.found),"\n")
        cat("Model, % Consistency Met=",c(this.m_label,p.consistency),"\n")
        }
       }
      }
  }
  
# Choice of subgroup  
out_hr <- out_Nsg <- out_Msg <- NULL
df_flag <- df_flag <- sg.harm <- sg.harm.id <- NULL

if(any.found>0){
result$Pcons <- as.numeric(result$Pcons)
result$hr <- as.numeric(result$hr)
result$N <- as.numeric(result$N)
result$E <- as.numeric(result$E)
result$K <- as.numeric(result$K)
result_new <- copy(result)

# sg_focus="hr"
# Option (details=TRUE) to display details and sg plot
sgdetails<-ifelse(plot.sg & sg_focus=="hr", TRUE, FALSE)
out_hr<-sg_consistency_out(df=df,result_new=result_new,sg_focus="hr",details=sgdetails,plot.sg=sgdetails,index.Z=index.Z,names.Z=names.Z,by.risk=by.risk,confs_labels=confs_labels)

# sg_focus="Nsg"
sgdetails<-ifelse(plot.sg & sg_focus %in% c("Nsg","Nsg_only"), TRUE, FALSE)
out_Nsg<-sg_consistency_out(df=df,result_new=result_new,sg_focus="Nsg",details=sgdetails,plot.sg=sgdetails,index.Z=index.Z,names.Z=names.Z,by.risk=by.risk,confs_labels=confs_labels)

# sg_focus="Msg"
sgdetails<-ifelse(plot.sg & sg_focus %in% c("Msg","Msg_only"), TRUE, FALSE)
out_Msg<-sg_consistency_out(df=df,result_new=result_new,sg_focus="Msg",details=sgdetails,plot.sg=sgdetails,index.Z=index.Z,names.Z=names.Z,by.risk=by.risk,confs_labels=confs_labels)

# Return "data" to duplicate last version
if(sg_focus=="hr"){
  df_flag <- out_hr$df_flag
  this.1 <- out_hr$sg.harm_label
  sg.harm.id <- out_hr$sg.harm.id
}
if(sg_focus %in% c("Nsg","Nsg_only")){
  df_flag <- out_Nsg$df_flag
  this.1 <- out_Nsg$sg.harm_label
  sg.harm.id <- out_Nsg$sg.harm.id
}
if(sg_focus %in% c("Msg","Msg_only")){
  df_flag <- out_Msg$df_flag
  this.1 <- out_Msg$sg.harm_label
  sg.harm.id <- out_Msg$sg.harm.id
}
if(details) cat("SG focus=",c(sg_focus),"\n")
}

if(details){
t.end<-proc.time()[3]
t.min<-(t.end-t.start)/60
cat("Subgroup Consistency Minutes=",c(t.min),"\n")
if(any.found>0){
cat("Subgroup found (FS)","\n")
}
if(any.found==0){
cat("NO subgroup found (FS)","\n")
  }
}
out<-list(out_hr=out_hr,out_Nsg=out_Nsg,out_Msg=out_Msg,df_flag=df_flag,sg.harm=this.1,sg.harm.id=sg.harm.id)
return(out)
}


# For plotting subgroups and the complement

plot.subgroup<-function(tte.name,event.name,treat.name,weight.name=NULL,
                        strata.name=NULL,sub1,sub1C,xloc1=NULL,xloc2=NULL,details=FALSE,show.logrank=FALSE,
                        ymin=0,cox.cex=0.6,prob.points=c(0,0.25,0.5,0.75,1.0),
                        cex_Yaxis=1,title1=NULL,title2=NULL,choose_ylim=TRUE,
                        exp.lab="Treat",con.lab="Control",legend.cex=0.60,risk.cex=0.60,
                        yloc1=0.6,yloc2=0.6,subid=NULL,byrisk=2,fix.rows=TRUE,show.med=FALSE,ylab="Survival"){
  
  
  if(!is.null(subid) & (!is.null(title1) | !is.null(title2))) stop("subid OR titles need to be specified (not both)")
  
  #if(is.null(weight.name)){
  #  sub1$wgt<-rep(1,nrow(sub1)) 
  #  sub1C$wgt<-rep(1,nrow(sub1C))
  #}
  
  #if(!is.null(wgt.name)){
  #  sub1$wgt<-sub1[,c(wgt.name)] 
  #  sub1C$wgt<-sub1C[,c(wgt.name)]
  #}
  
  # Set tpoints.add to span range of both subgroups
  aa<-max(sub1[,c(tte.name)])
  bb<-max(sub1C[,c(tte.name)])
  tpoints.add<-c(-1,max(aa,bb))
  
  if(fix.rows) par(mfrow=c(1,2))
  df<-sub1
  
  km.fit<-KM.plot.2sample.weighted(df=df, tte.name=tte.name,event.name=event.name, treat.name=treat.name,
                                   weight.name=weight.name,strata.name=strata.name,
                                   show_arm_legend=FALSE,
                                   risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
                                   cex_Yaxis=cex_Yaxis,choose_ylim=choose_ylim,
                                   stop.onerror=TRUE,Xlab="Months",Ylab=ylab,details=details,prob.points=prob.points,
                                   ymin=ymin,show.logrank=show.logrank,
                                   show.med=TRUE,show.cox=TRUE,cox.cex=cox.cex,med.cex=cox.cex,qlabel="med=")
  cpoints <- km.fit$cpoints
  
  m1<-round(km.fit$med.1,1)
  m0<-round(km.fit$med.2,1)
  
  if(is.null(xloc1)) xloc1<-c(diff(range(cpoints))/2)  
  
  #if(!show.med) legend(xloc1, yloc1, c(paste(c(exp.lab,eval(parse(text=m1))),collapse="; med="), paste(c(con.lab,eval(parse(text=m0))),collapse="; med=")), 
  #                     col = c("black","blue"), lwd = 2, bty = "n", cex=legend.cex)
  #if(show.med) legend(xloc1, yloc1, c(exp.lab, con.lab), col = c("black","blue"), lwd = 2, bty = "n")
  
  if(show.med){
    legend("top", c(paste(c(exp.lab,eval(parse(text=m1))),collapse="; med="), paste(c(con.lab,eval(parse(text=m0))),collapse="; med=")), 
           col = c("black","blue"), lwd = 2, bty = "n", cex=legend.cex)
  }
  if(!show.med) legend("top", c(exp.lab, con.lab), col = c("black","blue"), lwd = 2, bty = "n", cex=legend.cex)
  
  #legend("top",c(subid),bty="n",cex=0.7)
  if(!is.null(subid)) title(main=subid,cex.main=0.80)
  if(!is.null(title1)) title(main=title1,cex.main=0.80)
  
  df<-sub1C
  
  km.fit<-KM.plot.2sample.weighted(df=df, tte.name=tte.name,event.name=event.name, treat.name=treat.name,
                                   weight.name=weight.name,strata.name=strata.name,
                                   show_arm_legend=FALSE,
                                   risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
                                   cex_Yaxis=cex_Yaxis,choose_ylim=choose_ylim,
                                   stop.onerror=TRUE,Xlab="Months",Ylab=ylab,details=details,prob.points=prob.points,
                                   ymin=ymin,show.logrank=show.logrank,
                                   show.med=TRUE,show.cox=TRUE,cox.cex=cox.cex,med.cex=cox.cex,qlabel="med=")
  
  cpoints <- km.fit$cpoints
  m1<-round(km.fit$med.1,1)
  m0<-round(km.fit$med.2,1)
  
  if(is.null(xloc2)) xloc2<-c(diff(range(cpoints))/2)  
  
  if(show.med){
    legend("top", c(paste(c(exp.lab,eval(parse(text=m1))),collapse="; med="), paste(c(con.lab,eval(parse(text=m0))),collapse="; med=")), 
           col = c("black","blue"), lwd = 2, bty = "n",cex=legend.cex)
  }  
  if(!show.med) legend("top", c(exp.lab, con.lab), col = c("black","blue"), lwd = 2, bty = "n", cex=legend.cex)
  
  
  #if(!show.med) legend(xloc2, yloc2, c(paste(c(exp.lab,eval(parse(text=m1))),collapse="; med="), paste(c(con.lab,eval(parse(text=m0))),collapse="; med=")), 
  #                     col = c("black","blue"), lwd = 2, bty = "n",cex=legend.cex)
  #if(show.med) legend(xloc2, yloc2, c(exp.lab, con.lab), col = c("black","blue"), lwd = 2, bty = "n")
  
  
  if(!is.null(subid)) title(main="Complement",cex.main=0.80)
  if(!is.null(title2)) title(main=title2,cex.main=0.80)
}

