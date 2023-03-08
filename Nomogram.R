#install.packages("rms")
library(rms)
setwd("..")                        

riskFile="tcgaRisk.txt"
cliFile="tcgaClinical.txt"
outFile="tcga.Nomogram.pdf"
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
dd <- datadist(rt)
options(datadist="dd")
f <- cph(Surv(futime, fustat) ~ Age+riskScore, x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)
nom  <- nomgram(f, fun=list(function(x) surv(1, x), function(x) surv(5, x), function(x) surv(10, x)), 
    lp=F, funlabel=c("1-year survival", "5-year survival", "10-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.95, 0.9, 0.8, 0.7, 0.4, 0.2, 0.05))  
#nomfile=outFile,height=7.5,width=15)
plot(nom)
dev.off()


library(survival)
library(survminer)
library(timeROC)
library(rms)
setwd("..")    
bioROC=function(inputFile=null, rocFile=null, calFile=null){
  
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  
  ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
                 marker=rt$riskScore, cause=1,
                 weighting='aalen',
                 times=c(1,5,10), ROC=TRUE)
  pdf(file=rocFile,width=5,height=5)
  plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
  plot(ROC_rt,time=5,col='blue',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=10,col='red',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 10 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("green","blue","red"),lwd=2,bty = 'n')
  dev.off()
  
  
  pdf(file=calFile, width=5, height=5)
  
  f <- cph(Surv(futime, fustat) ~ riskScore, x=T, y=T, surv=T, data=rt, time.inc=1)
  cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
  plot(cal, xlim=c(0,1), ylim=c(0,1),
       xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)
  
  f <- cph(Surv(futime, fustat) ~ riskScore, x=T, y=T, surv=T, data=rt, time.inc=5)
  cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
  plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)
  
  f <- cph(Surv(futime, fustat) ~ riskScore, x=T, y=T, surv=T, data=rt, time.inc=10)
  cal <- calibrate(f, cmethod="KM", method="boot", u=10, m=(nrow(rt)/3), B=1000)
  plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
  legend('bottomright', c('1-year', '5-year', '10-year'),
         col=c("green","blue","red"), lwd=1.5, bty = 'n')
  dev.off()
}


bioROC(inputFile="risk.txt", rocFile=".ROC.pdf", calFile=".cal.pdf")

