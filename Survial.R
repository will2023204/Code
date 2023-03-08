#install.packages("survival")
#install.packages("survminer")
library(survival)
library("survminer")
setwd("..")            

bioSurvival=function(inputFile=null,outFile=null){
		diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
		pValue=1-pchisq(diff$chisq,df=1)
		pValue=signif(pValue,4)
		pValue=format(pValue, scientific = TRUE)
		fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
		
		surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=TRUE,
		           pval=paste0("p=",pValue),
		           pval.size=6,
		           risk.table=FALSE,
		           legend.labs=c("High risk", "Low risk"),
		           legend.title="Risk",
		           xlab="Time(years)",
		           break.time.by = 2,
		           risk.table.title="",
		           palette=c("red", "blue"),
		           risk.table.height=.25)
		pdf(file=outFile,onefile = FALSE,width = 7,height =6)
		print(surPlot)
		dev.off()
}
bioSurvival(inputFile="tcgaRisk.txt",outFile="tcgaRisk.pdf")
bioSurvival(inputFile="geoRisk.txt",outFile="geoRisk.pdf")


library(pheatmap)
bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null,heatmapFile=null){
  rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)  
  rt=rt[order(rt$riskScore),]                                           
  
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=riskScoreFile,width = 10,height = 3.5)
  plot(line, type="p", pch=20,
       ylab="Risk score",
       col=c(rep("green",lowLength),rep("red",highLength)) )
  abline(h=median(rt$riskScore),v=lowLength,lty=2)
  legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  dev.off()
  
  color=as.vector(rt$fustat)
  color[color==1]="red"
  color[color==0]="green"
  pdf(file=survStatFile,width = 10,height = 3.5)
  plot(rt$futime, pch=19,
       ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  
  rt1=rt[c(3:(ncol(rt)-2))]
  rt1=log2(rt1+1)
  rt1=t(rt1)
  annotation=data.frame(type=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  pdf(file=heatmapFile,width = 10,height = 4.5)
  pheatmap(rt1, 
           annotation=annotation, 
           cluster_cols = FALSE,
           fontsize_row=11,
           show_colnames = F,
           fontsize_col=3,
           color = colorRampPalette(c("green", "black", "red"))(50) )
  dev.off()
}
bioRiskPlot(inputFile="tcgaRisk.txt",riskScoreFile="tcga.riskScore.pdf",survStatFile="tcga.survStat.pdf",heatmapFile="tcga.heatmap.pdf")
