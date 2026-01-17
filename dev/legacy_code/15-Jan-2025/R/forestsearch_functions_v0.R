
# modified from ade4 package to handle single factor variable
acm.disjctif<-function (df) 
{
  acm.util.df <- function(i) {
    cl <- df[, i]
    cha <- names(df)[i]
    n <- length(cl)
    cl <- as.factor(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (unclass(cl) - 1)] <- 1
    dimnames(x) <- list(row.names(df), paste(cha, levels(cl), 
                                             sep = "."))
    return(x)
  }
  
  # For FAC(df) with only single variable
  acm.util.df2 <- function(i) {
    cl <- df[, i]
    cha <- colnames(df)[i]
    n <- length(cl)
    cl <- as.factor(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (unclass(cl) - 1)] <- 1
    dimnames(x) <- list(row.names(df), paste(cha, levels(cl), 
                                             sep = "."))
    return(x)    
  }
  
  if(!is.null(ncol(df)) & ncol(df)>1) G <- lapply(1:ncol(df), acm.util.df)
    
  if(!is.null(ncol(df)) & ncol(df)==1) G <- lapply(1, acm.util.df2)
    
  G <- data.frame(G, check.names = FALSE)
  return(G)
}


dummy2 <- function(df) {  
  
  NUM <- function(dataframe)dataframe[,sapply(dataframe,is.numeric)]
  FAC <- function(dataframe)dataframe[,sapply(dataframe,is.factor)]
  
  if (is.null(ncol(NUM(df)))){
    DF <- data.frame(NUM(df), acm.disjctif(FAC(df)))
    names(DF)[1] <- colnames(df)[which(sapply(df, is.numeric))]
  } else {
    if (!is.null(ncol(FAC(df))) & ncol(FAC(df))>0) DF <- data.frame(NUM(df), acm.disjctif(FAC(df)))
    if (!is.null(ncol(FAC(df))) | ncol(FAC(df))==0) DF <- data.frame(NUM(df))
    
    if (is.null(ncol(FAC(df)))){
      temp<-as.matrix(FAC(df))  
      colnames(temp)[1]<-colnames(df)[which(sapply(df, is.factor))]
      df.fac<-acm.disjctif(temp)  
      DF <- data.frame(NUM(df), df.fac)  
    }
  }
  return(DF)
} 


dummy <- function(df) {  
  
  NUM <- function(dataframe)dataframe[,sapply(dataframe,is.numeric)]
  FAC <- function(dataframe)dataframe[,sapply(dataframe,is.factor)]
  
  if (is.null(ncol(NUM(df)))){
    DF <- data.frame(NUM(df), acm.disjctif(FAC(df)))
    names(DF)[1] <- colnames(df)[which(sapply(df, is.numeric))]
  } else {
    if (!is.null(ncol(FAC(df)))) DF <- data.frame(NUM(df), acm.disjctif(FAC(df)))
     if (is.null(ncol(FAC(df)))){
      temp<-as.matrix(FAC(df))  
      colnames(temp)[1]<-colnames(df)[which(sapply(df, is.factor))]
      df.fac<-acm.disjctif(temp)  
      DF <- data.frame(NUM(df), df.fac)  
    }
  }
  return(DF)
} 





ztrail<-function(kk){
  ii<-1
  zz<-kk
  while(zz%%2==0){
    ii<-ii+1
    zz<-zz/2
  }
  return(ii)
}

one.zero<-function(x){
  if (x==1){
    x<-0}
  else {
    x<-1}
  return(x)
}

# Consistency splitting

split.df<-function(df,r.train){
n.split<-round(r.train*nrow(df),0)
data.id<-c(1:nrow(df))
id.split1<-sample(data.id,size=n.split,replace=F)
df.split1<-df[id.split1,]
id.split2<-setdiff(data.id,id.split1)
df.split2<-df[id.split2,]
df.train<-df.split1
df.test<-df.split2
if(nrow(df.train)+nrow(df.test)!=nrow(df)) stop("Problem with cross-validation splitting (n(train)+n(est)!=n(all))")
return(list(df.train=df.train,df.test=df.test))
}

# The following is not used but keeping just-in-case

# For categorizing continuous variables by quantiles
cutq<-function(x,labels=c("q1","q2","q3","q4")){
  q.25<-quantile(x,0.25)
  q.75<-quantile(x,0.75)
  if(min(x)!=q.25 & max(x)!=q.75){
    cut.x<-cut(x,quantile(x,c(0,0.25,0.5,0.75,1)),include.lowest=TRUE,labels=labels)
  }
  
  if(min(x)==q.25 & max(x)!=q.75){
    cut.x<-cut(x,quantile(x,c(0.25,0.5,0.75,1)),include.lowest=TRUE,labels=c("2","3","4"))
  }
  
  if(min(x)==q.25 & max(x)==q.75){
    cut.x<-cut(x,quantile(x,c(0.25,0.5,1)),include.lowest=TRUE,labels=c("2","3"))
  }
  
  if(min(x)!=q.25 & max(x)==q.75){
    cut.x<-cut(x,quantile(x,c(0,0.25,0.5,1)),include.lowest=TRUE,labels=c("1","2","3"))
  }
  
  if(length(unique(cut.x))==2) cut.x<-as.numeric(cut.x) 
  return(cut.x)
}

cutmed<-function(x,direction=">med"){
  if(direction==">med"){
    cut.x<-ifelse(x>median(x),1,0)
  }
  if(direction=="<=med"){
    cut.x<-ifelse(x<=median(x),1,0)
  }
  return(cut.x)
}

cutmed.factor<-function(x,direction=">med",labels=c("1","2")){
  if(direction==">med"){
    cut.x<-factor(ifelse(x>median(x),1,0),levels=0:1,labels=labels)
  }
  if(direction=="<=med"){
    cut.x<-factor(ifelse(x<=median(x),1,0),levels=0:1,labels=labels)
  }
  return(cut.x)
}

cut1.factor<-function(x,cutpoint=median(x),direction="GT",labels=c("1","2")){
  if(direction=="GT"){
    cut.x<-factor(ifelse(x>cutpoint,1,0),levels=0:1,labels=labels)
  }
  if(direction=="LTE"){
    cut.x<-factor(ifelse(x<=cutpoint,1,0),levels=0:1,labels=labels)
  }
  return(cut.x)
}

cut2.factor<-function(x,cutpoints=c(quantile(x,c(0,0.25,0.5,1))),labels=c("1","2","3")){
  cut.x<-cut(x,cutpoints,include.lowest=TRUE,labels=labels)
  return(cut.x)
}


cut1.num<-function(x,cutpoint=median(x),direction="LTE"){
  if(direction=="GT"){
    cut.x<-ifelse(x>cutpoint,1,0)
  }
  if(direction=="LTE"){
    cut.x<-ifelse(x<=cutpoint,1,0)
  }
  return(cut.x)
}


# cut.factor and cut.num previously in sim_aft_gbsg function files
cut.factor<-function(x,cutpoint=median(x),direction="LTE",labels=c("1","2")){
  if(direction=="GT"){
    cut.x<-factor(ifelse(x>cutpoint,1,0),levels=0:1,labels=labels)
  }
  if(direction=="LTE"){
    cut.x<-factor(ifelse(x<=cutpoint,1,0),levels=0:1,labels=labels)
  }
  return(cut.x)
}

cut.num<-function(x,cutpoint=median(x),direction="LTE"){
  if(direction=="GT"){
    cut.x<-ifelse(x>cutpoint,1,0)
  }
  if(direction=="LTE"){
    cut.x<-ifelse(x<=cutpoint,1,0)
  }
  return(cut.x)
}

