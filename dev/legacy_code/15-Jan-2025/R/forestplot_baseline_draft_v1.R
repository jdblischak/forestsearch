require(forestploter)
forestplot_baseline <- function(df,confounders.name,arrow_text=c("favors Treatment","Control"),E.name=c("Treat"),C.name=c("Control"),outcome.name,event.name,treat.name,footnote_text=NULL,
xlim = c(0.25, 2.0), ticks_at = c(0.65, 1.0, 1.3),title_text=NULL){
                            
title_text <- NULL

  # ITT 
  res_itt <- SG_HRtable2(dfa=df,df.OOB=NULL,fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,
                         sg_name="ITT",E.name=E.name,C.name=C.name)
  # Add separator header 
  zz <- res_itt  
  zz$Subgroup <- "                 "
  zz[,c(E.name,C.name)] <- ""
  zz[,c("est","low","hi","se")] <- NA
  
  
  dt <- rbind(res_itt,zz)

  for(sg_index in 1:length(confounders.name)){
    
  sg_name <- confounders.name[sg_index]
  
  cov <- df[,sg_name]
  cov_values <- sort(unique(cov),decreasing=TRUE)
 
 if(length(cov_values)>=3){
  cut <- c(quantile(cov,c(0.5)))

  thiscut <- paste0(sg_name,"<=",cut)  
  dfsg <- subset(df,eval(parse(text=thiscut)))

  resA <- SG_HRtable2(dfa=dfsg,df.OOB=NULL, sg_name=thiscut,
  fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,E.name=E.name,C.name=C.name)
  
  
  thiscut <- paste0(sg_name,">",cut)  
  dfsg <- subset(df,eval(parse(text=thiscut)))
  
  resB <- SG_HRtable2(dfa=dfsg,df.OOB=NULL, sg_name=thiscut,
                      fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,E.name=E.name,C.name=C.name)
  dt <- rbind(dt,resA,resB,zz)
  rm("resA","resB")
 }
  
  if(length(cov_values)==2){
  
    thiscut <- paste0(sg_name,"==",cov_values[1])  
    dfsg <- subset(df,eval(parse(text=thiscut)))
    
    resA <- SG_HRtable2(dfa=dfsg,df.OOB=NULL, sg_name=thiscut,
                       fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,E.name=E.name,C.name=C.name)
    
    
    thiscut <- paste0(sg_name,"==",cov_values[2])  
    dfsg <- subset(df,eval(parse(text=thiscut)))
    
    resB <- SG_HRtable2(dfa=dfsg,df.OOB=NULL, sg_name=thiscut,
                       fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,E.name=E.name,C.name=C.name)
    dt <- rbind(dt,resA,resB,zz)
    rm("resA","resB")
  }
  
  }  
  
  dt <- dt[-nrow(dt),]
  
  sg_colors <- c("white")
  
  tm <- forest_theme(core=list(fg_params=list(hjust = 1, x = 0.9),
                               bg_params=list(fill = sg_colors)),
                     colhead=list(fg_params=list(hjust=0.5, x=0.5)),
                     footnote_gp = gpar(cex = 0.7, fontface = "italic", col = "blue"))
  
  dt$` ` <- paste(rep(" ", 20), collapse = " ")
  
  # Create a confidence interval column to display
  dt$`HR (95% CI)` <- ifelse(is.na(dt$se), "",
                             sprintf("%.2f (%.2f to %.2f)", dt$est, dt$low, dt$hi))
  
  
  p <- forest(dt[,c(1:3, 8:9)],
              title=title_text,
              est = as.numeric(dt$est),
              lower = as.numeric(dt$low), 
              upper = as.numeric(dt$hi),
              sizes = as.numeric(dt$se),
              ci_column = 4,
              ref_line = 1,
              arrow_lab = arrow_text,
              xlim = xlim,
              ticks_at = ticks_at,
              footnote = footnote_text, theme=tm)
  
  
return(list(dt=dt,p=p))
}
  
 