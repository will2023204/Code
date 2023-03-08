
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("car", "ridge", "preprocessCore", "genefilter", "sva"))

#install.packages("ggpubr")

library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
expFile="symbol.txt"     
riskFile="risk.txt"      
drug="Cisplatin"        
setwd("..")  

rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=t(data)

senstivity=pRRopheticPredict(data, drug, selection=1)
senstivity=senstivity[senstivity!="NaN"]
senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(risk), names(senstivity))
risk=risk[sameSample, "risk",drop=F]
senstivity=senstivity[sameSample]
rt=cbind(risk, senstivity)

rt$risk=factor(rt$risk, levels=c("low", "high"))
type=levels(factor(rt[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

boxplot=ggboxplot(rt, x="risk", y="senstivity", fill="risk",
			      xlab="Risk",
			      ylab=paste0(drug, " senstivity (IC50)"),
			      legend.title="Risk",
			      palette=c("#4DBBD5", "#E64B35")
			     )+ 
	stat_compare_means(comparisons=my_comparisons)
pdf(file=paste0(drug, ".pdf"), width=5, height=4.5)
print(boxplot)
dev.off()


library(limma)
library(ggpubr)
riskFile="cohort risk.txt"      
cliFile="cohort response.txt"           
setwd("..")     

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$riskScore[risk$riskScore>quantile(risk$riskScore,0.99)]=quantile(risk$riskScore,0.99)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

samSample=intersect(row.names(risk), row.names(cli))
risk=risk[samSample,"riskScore",drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk, cli)
rt=rt[,-3]
for(clinical in colnames(rt)[2:ncol(rt)]){
  data=rt[c("riskScore", clinical)]
  colnames(data)=c("riskScore", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  group=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  boxplot=ggboxplot(data, x="clinical", y="riskScore", color="clinical",
                    xlab="",
                    ylab="Risk score",
                    legend.title=clinical,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  pdf(file=paste0("cohort", clinical, ".pdf"), width=6, height=5)
  print(boxplot)
  dev.off()
}

library(plyr)
library(ggplot2)
library(ggpubr)
scoreFile="cohort risk.txt"      
cliFile="cohort response.txt"     
trait="Response"                 
setwd("..")   
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(score), row.names(cli))
rt=cbind(score[sameSample,], cli[sameSample,])
bioCol=c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85")
colnames(rt)[30]="Response"
trait=colnames(rt[ncol(rt)])
bioCol=bioCol[1:length(unique(rt[,trait]))]

rt1=rt[,c(trait, "risk")]
colnames(rt1)=c("trait", "risk")
df=as.data.frame(table(rt1))
df=ddply(df, .(risk), transform, percent = Freq/sum(Freq) * 100)
df=ddply(df, .(risk), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
p=ggplot(df, aes(x = factor(risk), y = percent, fill = trait)) +
  geom_bar(position = position_stack(), stat = "identity", width = .6) +
  scale_fill_manual(values=bioCol)+
  xlab("Risk group")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
  #coord_flip()+
  theme_bw()
pdf(file=paste0(trait,".pdf"), width=3.5, height=4)
print(p)
dev.off()

library(limma)
expFile="symbol.txt"         
riskFile="tcgaRisk.txt"       
logFCfilter=1                
fdrFilter=0.05               
setwd("..")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=t(data)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(colnames(data), row.names(risk))
data=data[,sameSample]
risk=risk[sameSample,]

riskLow=risk[risk$risk=="low",]
riskHigh=risk[risk$risk=="high",]
dataLow=data[,row.names(riskLow)]
dataHigh=data[,row.names(riskHigh)]
data=cbind(dataLow,dataHigh)
data=data[rowMeans(data)>1,]
conNum=ncol(dataLow)
treatNum=ncol(dataHigh)
Type=c(rep(1,conNum), rep(2,treatNum))

outTab=data.frame()
for(i in row.names(data)){
  rt=data.frame(expression=data[i,], Type=Type)
  wilcoxTest=wilcox.test(expression ~ Type, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,lowMean=conGeneMeans,highMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff, file="riskDiff.txt", sep="\t", row.names=F, quote=F)
