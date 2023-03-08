#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("scales")
#install.packages("ggplot2")
#install.packages("ggtext")

library(limma)
library(scales)
library(ggplot2)
library(ggtext)
riskFile="tcgaRisk.txt"     
immFile="infiltration.csv"     
setwd("..")     

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

immune=read.csv(immFile, header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))
immune=avereps(immune)

sameSample=intersect(row.names(risk), row.names(immune))
risk=risk[sameSample, "riskScore"]
immune=immune[sameSample,]

x=as.numeric(risk)
outTab=data.frame()
for(i in colnames(immune)){
	y=as.numeric(immune[,i])
	corT=cor.test(x, y, method="spearman")
	cor=corT$estimate
	pvalue=corT$p.value
	if(pvalue<0.05){
		outTab=rbind(outTab,cbind(immune=i, cor, pvalue))
	}
}
write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)

corResult=read.table("corResult.txt", head=T, sep="\t")
corResult$Software=sapply(strsplit(corResult[,1],"_"), '[', 2)
corResult$Software=factor(corResult$Software,level=as.character(unique(corResult$Software[rev(order(as.character(corResult$Software)))])))
b=corResult[order(corResult$Software),]
b$immune=factor(b$immune,levels=rev(as.character(b$immune)))
colslabels=rep(hue_pal()(length(levels(b$Software))),table(b$Software))    
pdf(file="cor.pdf", width=8, height=8)       #????ͼƬ
ggplot(data=b, aes(x=cor,y=immune, color=Software))+
	labs(x="Correlation coefficient",y="Immune cell")+
	geom_point(size=4.1)+
	theme(panel.background=element_rect(fill="white",size=1,color="black"),
	      panel.grid=element_line(color="grey75",size=0.5),
	      axis.ticks = element_line(size=0.5),
	      axis.text.y = ggtext::element_markdown(colour=rev(colslabels)))
dev.off()



#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library("limma")         
expFile="symbol.txt"    
setwd("..")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

v=voom(data, plot=F, save.plot=F)
out=v$E
out=rbind(ID=colnames(out), out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)   

source("13.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000, QN=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")

library(limma)
library(reshape2)
library(ggpubr)

riskFile="tcgaRisk.txt"            
immFile="CIBERSORT-Results.txt"    
pFilter=0.05                     
setwd("..")    
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])

group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
rt=cbind(data[sameSample,,drop=F], risk[sameSample,"Risk",drop=F])
rt=rt[order(rt$Risk, decreasing=T),]
conNum=nrow(rt[rt$Risk=="low",])
treatNum=nrow(rt[rt$Risk=="high",])

data=t(rt[,-ncol(rt)])
pdf("barplot.pdf", height=20, width=36)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1], ybottom = -0.01, xright = a1[conNum], ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"Low risk",cex=2)
rect(xleft = a1[conNum], ybottom = -0.01, xright =a1[length(a1)] , ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"High risk",cex=2)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()


data=rt
data=melt(data, id.vars=c("Risk"))
colnames(data)=c("Risk", "Immune", "Expression")
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("low","high"))
bioCol=c("#4DBBD5","#E64B35","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Risk",
                  xlab="",
                  ylab="Fraction",
                  legend.title="Risk",
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()


immFile="CIBERSORT-Results.txt"     
riskFile="tcgaRisk.txt"            
setwd("..")     

immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])

group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,3:(ncol(risk)-2),drop=F]

outTab=data.frame()
for(immune in colnames(data)){
  for(gene in colnames(risk)){
    x=as.numeric(data[,immune])
    y=as.numeric(risk[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, Immune=immune, cor, text, pvalue))
  }
}

outTab$cor=as.numeric(outTab$cor)
pdf(file="cor.pdf", width=8, height=6)
ggplot(outTab, aes(Gene, Immune)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() +    
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),   
        axis.text.y = element_text(size = 10, face = "bold")) +       #y??????
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   
  scale_x_discrete(position = "bottom")     
dev.off()


library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)
library(estimate)
inputFile="symbol.txt"      
setwd("..")       

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

out=rbind(ID=colnames(data),data)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")

scores=read.table("estimateScore.gct", skip=2, header=T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.", "\\-", rownames(scores))
out=rbind(ID=colnames(scores), scores)
write.table(out, file="TMEscores.txt", sep="\t", quote=F, col.names=F)


riskFile="tcgaRisk.txt"          
estimateFile="TMEscores.txt"      
setwd("..")     

Risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
Risk$risk=factor(Risk$risk, levels=c("low","high"))

score=read.table(estimateFile, header=T, sep="\t", check.names=F, row.names=1)
score=score[,1:3]
rownames(score)=gsub("(.*?)\\_(.*?)", "\\2", rownames(score))
score=score[row.names(Risk),,drop=F]

rt=cbind(Risk[,"risk",drop=F], score)

data=melt(rt, id.vars=c("risk"))
colnames(data)=c("Risk", "scoreType", "Score")
write.table(data,".txt",quote = F,sep = "\t",col.names = T,row.names = F)
p=ggviolin(data, x="scoreType", y="Score", fill = "Risk",
           xlab="",
           ylab="TME score",
           legend.title="Risk",
           add = "boxplot", add.params = list(color="white"),
           palette = c("blue","red"), width=1)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Risk),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

pdf(file="vioplot.pdf", width=6, height=5)
print(p1)
dev.off()


library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
riskFile="tcgaRisk.txt"       
RNAssFile="StemnessScore.tsv"          
setwd("..")     

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

RNAss=read.table(RNAssFile, header=T, sep="\t",check.names=F, row.names=1)
RNAss=t(RNAss[1,,drop=F])
rownames(RNAss)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(RNAss))
RNAss=avereps(RNAss)

sameSample=intersect(row.names(risk), row.names(RNAss))
risk=risk[sameSample,c("riskScore","risk"),drop=F]
RNAss=RNAss[sameSample,,drop=F]
data=cbind(RNAss, risk)
write.table(data,".txt",quote = F,sep = "\t",col.names = T,row.names = F)
xlab="riskScore"
ylab="RNAss"
outFile="RNAss.cor.pdf"
x=as.numeric(data[,xlab])
x[x>quantile(x,0.99)]=quantile(x,0.99)
y=as.numeric(data[,ylab])
df1=as.data.frame(cbind(x,y))
p1=ggplot(df1, aes(x, y)) + 
  xlab("Risk score") + ylab(ylab)+ ylim(0,0.7)+
  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
pdf(file=outFile, width=5.2, height=5)
print(p2)
dev.off()


