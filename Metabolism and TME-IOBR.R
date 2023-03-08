library(IOBR)
library(EPIC)
library(estimate) 
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(survival)
library(UCSCXenaTools)
setwd("..")

eset_stad<-read.table("BRCA-TCGA-TPM.txt",sep = "\t",header = T,row.names = 1,check.names = F)

eset_stad[1:5, 1:5]
names(signature_metabolism)[1:20]
eset_stad[1:5, 1:5]
nrow(eset_stad)
ncol(eset_stad)
#calculate the signature gene sets derived from GO, KEGG, HALLMARK and REACTOME datasets
sig_hallmark<-calculate_sig_score(pdata           = NULL,
                                  eset            = eset_stad,
                                  signature       = hallmark,
                                  method          = "ssgsea",
                                  mini_gene_count = 2)

sig_hallmark[1:5,1:10]
write.table(sig_hallmark,"BRCA-hallmark_signature.txt",quote = F,sep = "\t",row.names = F,col.names = T)
#Calculate TME associated signatures-(through PCA method).
sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset_stad,
                             signature       = signature_tme,
                             method          = "pca",
                             mini_gene_count = 2)
sig_tme[1:5, 1:10]

#Estimate TME associated signatures-(through ssGSEA method).
sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset_stad,
                             signature       = signature_tme,
                             method          = "ssgsea",
                             mini_gene_count = 5)
sig_tme[1:5, 1:10]
write.table(sig_tme,"BRCA-TME_signature.txt",quote = F,sep = "\t",row.names = T,col.names = T)

#Evaluate metabolism related signatures.
sig_meta<-calculate_sig_score(pdata           = NULL,
                              eset            = eset_stad,
                              signature       = signature_metabolism,
                              method          = "ssgsea",
                              mini_gene_count = 2)
sig_meta[1:5, 1:10]
write.table(sig_meta,"BRCA-Metab_signature.txt",quote = F,sep = "\t",row.names = F,col.names = T)

#Evaluate Signatures associated with biomedical basic research: such as m6A and exosomes.
sig_basic<-calculate_sig_score(pdata           = NULL,
                              eset            = eset_stad,
                              signature       = signature_tumor,
                              method          = "ssgsea",
                              mini_gene_count = 2)
sig_basic[1:5, 1:10]
write.table(sig_basic,"BRCA-Basic_signature.txt",quote = F,sep = "\t",row.names = F,col.names = T)

# Analyze all collected signature scores (integrating three methods: PCA, ssGSEA and z-score).
sig_res<-calculate_sig_score(pdata           = NULL,
                             eset            = eset_stad,
                             signature       = signature_collection,
                             method          = "integration",
                             mini_gene_count = 2)
write.table(sig_res,"Research_signature_intergration.txt",quote = F,sep = "\t",row.names = F,col.names = T)
sig_res_gsea<-calculate_sig_score(pdata           = NULL,
                             eset            = eset_stad,
                             signature       = signature_collection,
                             method          = "ssgsea",
                             mini_gene_count = 2)
write.table(sig_res_gsea,"Research_signature_ssgsea.txt",quote = F,sep = "\t",row.names = F,col.names = T)
# signature group
sig_group$tumor_signature
a=as.matrix(a)
write.table(a,"signature_group.txt",quote = F,sep = "\t",row.names = T,col.names = T)

#citation of signatures
signature_collection_citation[1:20,]
a=signature_collection_citation
write.table(a,"signature_citation.txt",quote = F,sep = "\t",row.names = F,col.names = T)



library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)
library(gplots)
rm(list=ls())
expFile="symbol.txt"      
clusterFile="tcgaRisk.txt"      
setwd("..")      #???ù???Ŀ¼

gsvaCluster=read.table("TME_Basic_hallmark.txt",sep = "\t",header = T,check.names = F,row.names = 1)
Project=c(rep("TCGA",1023))
gsvaCluster=cbind(gsvaCluster, Project)
adj.P.Val.Filter=0.05
allType=as.vector(gsvaCluster$risk)
comp=combn(levels(factor(allType)), 2)
for(i in 1:ncol(comp)){
  con=gsvaCluster[gsvaCluster$risk==comp[2,i],]
  treat=gsvaCluster[gsvaCluster$risk==comp[1,i],]
  data=rbind(con, treat)
  Type=as.vector(data$risk)
  ann=data[,c(ncol(data), (ncol(data)-1))]
  data=t(data[,-c((ncol(data)-1), ncol(data))])
  design=model.matrix(~0+factor(Type))
  colnames(design)=levels(factor(Type))
  fit=lmFit(data, design)
  contrast=paste0(comp[1,i], "-", comp[2,i])
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)

  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".TME_basic_hallmark-all.txt"), sep="\t", quote=F, col.names=F)
  
  diffSig=allDiff[with(allDiff, adj.P.Val < adj.P.Val.Filter ), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".TME_basic_hallmark-all-diff.txt"), sep="\t", quote=F, col.names=F)
  
  bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  ann_colors=list()
  m6aCluCol=bioCol[1:length(levels(factor(allType)))]
  names(m6aCluCol)=levels(factor(allType))
  ann_colors[["risk"]]=m6aCluCol[c(comp[1,i], comp[2,i])]
  
  termNum=..    
  will=intersect(rownames(data),rownames(diffSig))
  diffSig=diffSig[will,]
  diffTermName=as.vector(rownames(diffSig))
  diffLength=length(diffTermName)
  if(diffLength<termNum){termNum=diffLength}
  hmGene=diffTermName[1:termNum]
  hmExp=data[hmGene,]
  pdf(file=paste0(contrast,".EME_Basic_Hallmark-heatmap.pdf"), width=9.6, height=6)
  pheatmap(hmExp, 
           annotation=ann,
           annotation_colors = ann_colors,
           color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
           #col = bluered(256),
           cluster_cols =F,
           cluster_rows = F,
           show_colnames = F,
           #gaps_col=as.vector(cumsum(table(Type))),
           gaps_row = as.vector(c(6,22)),
           scale="row",
           fontsize = 10,
           fontsize_row=10,
           fontsize_col=10)
  dev.off()
}


library(limma)
library(corrplot)
library(ggpubr)
library(ggExtra)
expFile=".txt"   
riskFile="tcgaRisk.txt"       
setwd("..")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data),row.names(risk))
rt1=cbind(data[sameSample,,drop=F],risk[sameSample,"riskScore",drop=F])

M=cor(rt1)
res1=cor.mtest(rt1, conf.level = 0.95)

pdf(file="COR.pdf", width=8, height=8)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=0.8, pch=T,
         #p.mat = res1$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()


