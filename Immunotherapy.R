
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
expFile="symbol.txt"      
riskFile="tcgaRisk.txt"       
geneFile="checkpoint.txt"      
setwd("..")    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
	
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])
data=log2(data+1)

group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(data))
data=avereps(data)
	
risk=read.table(riskFile, sep="\t", header=T, check.names=F, row.names=1)
sameSample=intersect(row.names(data),row.names(risk))
rt1=cbind(data[sameSample,],risk[sameSample,])
rt1=rt1[,c(sameGene,"risk")]
sigGene=c()
for(i in colnames(rt1)[1:(ncol(rt1)-1)]){
	if(sd(rt1[,i])<0.001){next}
	wilcoxTest=wilcox.test(rt1[,i] ~ rt1[,"risk"])
	pvalue=wilcoxTest$p.value
	if(wilcoxTest$p.value<0.05){
		sigGene=c(sigGene, i)
	}
}
sigGene=c(sigGene, "risk")
rt1=rt1[,sigGene]

rt1=melt(rt1,id.vars=c("risk"))
colnames(rt1)=c("risk","Gene","Expression")
group=levels(factor(rt1$risk))
rt1$risk=factor(rt1$risk, levels=c("low","high"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}

boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill="risk",
				  xlab="",
				  ylab="Gene expression",
				  legend.title="Risk",
				  width=0.8,
				  palette = c("#4DBBD5", "#E64B35") )+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=risk),
	method="wilcox.test",
	symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
	
pdf(file="checkpoint.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()


library(limma)
library(ggpubr)
tideFile="tide.txt"          
riskFile="tcgaRisk.txt"   
setwd("..")   
tide=read.table(tideFile, header=T, sep="\t", check.names=F, row.names=1)
group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tide=tide[group==0,]
tide=as.matrix(tide)
row.names(tide)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
tide=avereps(tide)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "Risk", drop=F]
data=cbind(tide, risk)

data$Risk=ifelse(data$Risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("Low-risk", "High-risk"))
group=levels(factor(data$Risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

for(i in colnames(data)[1:(ncol(data)-1)]){
  gg1=ggviolin(data, x="Risk", y=i, fill = "Risk", 
               xlab="", ylab=i,
               palette=c("#0066FF","#FF0000"),
               legend.title="Risk",
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  
  pdf(file=paste0("violin.", i, ".pdf"), width=6, height=5)
  print(gg1)
  dev.off()
}


library(reshape2)
library(ggpubr)
riskFile="tcgaRisk.txt"          
estimateFile="IPS score.txt"     
setwd("..")   

Risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
Risk$risk=factor(Risk$risk, levels=c("low","high"))

score=read.table(estimateFile, header=T, sep="\t", check.names=F, row.names=1)
score=score[,1:4]
samesample=intersect(rownames(Risk),rownames(score))
Risk=Risk[samesample,]
score=score[samesample,]
#rownames(score)=gsub("(.*?)\\_(.*?)", "\\2", rownames(score))
score=score[row.names(Risk),,drop=F]
rt=cbind(score,Risk[,"risk",drop=F])
write.table(rt,"risk-ips.txt",quote = F,sep = "\t",row.names = T,col.names = T)
data=melt(rt, id.vars=c("risk"))
colnames(data)=c("Risk", "scoreType", "Score")

p=ggviolin(data, x="scoreType", y="Score", fill = "Risk",
           xlab="",
           ylab="IPS score",
           legend.title="Risk",
           add = "boxplot", add.params = list(color="black"),
           palette = c("#4DBBD5","#E64B35"), width=1)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Risk),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

pdf(file="vioplot.pdf", width=6, height=5)
print(p1)
dev.off()


library(limma)
library(GSVA)
library(GSEABase)
library(ggpubr)
library(reshape2)

expFile="symbol.txt"          
gmtFile="immune.gmt"           
riskFile="tcgaRisk.txt"      
setwd("..")      

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]

geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file="immFunScore.txt", sep="\t", quote=F, col.names=F)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,c("Risk","riskScore"),drop=F]
rt1=cbind(data, risk)
write.table(data,"riskscore-immunefunction.txt",quote = F,sep = "\t",col.names = T,row.names = T)
data=melt(rt1, id.vars=c("Risk"))
colnames(data)=c("Risk","Type","Score")
data$Risk=factor(data$Risk, levels=c("low","high"))
p=ggboxplot(data, x="Type", y="Score", color = "Risk",
            xlab="",ylab="Score",add = "none",palette = c("#4DBBD5","#E64B35") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")

pdf(file="immFunction.pdf", width=8, height=6)
print(p)
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