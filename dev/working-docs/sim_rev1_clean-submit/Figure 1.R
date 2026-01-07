
rm(list=ls())
codepath<-c("/media/larryleon/My Projects/GitHub/Forest-Search/R/")
source(paste0(codepath,"source_forestsearch_v0.R"))
source_fs_functions(file_loc=codepath)

library(survival)
# This file is created in applications output
# Transfered from output to paper directory
load(file="../applications/simulated_ex1/output/sim-FS4_N1k_Noise=3_hrH=2_sim99.Rdata")
# Analysis dataset
df <- fs.est$df.est

# Use Rstudio file conversion button
#tiff("Figure 1.tiff", width=480, height=480, compression="lzw", res=600) 

tte.name<-c("y.sim")
event.name<-c("event.sim")
treat.name<-c("treat")

byrisk <- 6

risk.cex<-0.95
legend.cex<-0.90

# Add all events
evs<-sort(unique(df$y.sim[which(df$event.sim==1)]))
tpoints.add<-c(-1,evs,max(df$y.sim))

# gray40 also good for grey

kmH.fit<-KM.plot.2sample.weighted(Y=df[,c(tte.name)],E=df[,c(event.name)],Treat=df[,c(treat.name)],
                                  risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
                                  risk_offset=0.125, risk_delta=0.05,censor.cex=0.8,
                                  choose_ylim=FALSE,
                                  stop.onerror=TRUE, check.KM=TRUE,
                                  Xlab="Months",Ylab="Recurrence survival",details=FALSE,
                                  show.ticks=TRUE,
                                  col.1="black", col.2="blue",
                                  ltys=c(1,1),lwds=c(1,1),
                                  plotit=TRUE,
                                  show.logrank=FALSE,show.med=FALSE,show.cox=FALSE)
legend("top",legend=c("Treatment","Control"),col=c("black","blue"), lty=c(1,1),lwd=1,bty="n")

#dev.off()

