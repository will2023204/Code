

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

expFile="symbol.txt"          
riskFile="tcgaRisk.txt"     
gmtFile="c5.go.v7.4.symbols.gmt"     
setwd("..")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))

Risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
data=data[,row.names(Risk)]

dataL=data[,row.names(Risk[Risk[,"risk"]=="low",])]
dataH=data[,row.names(Risk[Risk[,"risk"]=="high",])]
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
logFC=sort(logFC,decreasing=T)
genes=names(logFC)
a=cbind(genes,logFC)
write.table(a,"risk-diff.txt",quote = F,sep = "\t",row.names = F,col.names = T)
gmt=read.gmt(gmtFile)

kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
write.table(kkTab,file="GSEA-GO.result.txt",sep="\t",quote=F,row.names = F)
	
termNum=5   
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
	showTerm=row.names(kkUp)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in high risk group")
	pdf(file="GSEA-GO.highRisk.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}

termNum=5    
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in low risk group")
	pdf(file="GSEA-GO.lowRisk.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}


