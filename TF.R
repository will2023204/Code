setwd("..")                 
fdrFilter=0.01  
logFCfilter=0.585

rt=read.table("TCGA.all.txt",sep="\t",header=T,check.names=F,row.names=1)
diffExp=read.table("TCGA.DiffExp.txt",sep="\t",header=T,check.names=F,row.names=1)
gene=read.table("TF.txt",sep="\t",header=F)
TF_DiffAll=rt[intersect(gene[,1],row.names(rt)),]
TF_DiffGene=intersect(gene[,1],row.names(diffExp))
hmExp=diffExp[TF_DiffGene,]
TF_DiffResult=TF_DiffAll[TF_DiffGene,]

write.table(TF_DiffAll,file="TF-ALL.xls",sep="\t",col.names=F,quote=F)
TF_DiffResult=rbind(ID=colnames(TF_DiffResult),TF_DiffResult)
write.table(TF_DiffResult,file="TFdiff.xls",sep="\t",col.names=F,quote=F)

eneExp=rbind(ID=colnames(hmExp),hmExp)
write.table(TF_GeneExp,file="TFgeneExp.txt",sep="\t",col.names=F,quote=F)

corFilter=0.3              
pvalueFilter=0.001         
setwd("..")                  

TF = read.table("TFgeneExp.txt", row.names=1 ,header=T,sep="\t",check.names=F)  
immuneGene = read.table("tcgaUniSigExp.txt", row.names=1 ,header=T,sep="\t",check.names=F)   
immuneGene=t(immuneGene[,3:ncol(immuneGene)])
rownames(immuneGene)=gsub("\\|","\\-",rownames(immuneGene))

colnames(TF)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(TF))
sameSample=intersect(colnames(TF),colnames(immuneGene))
TF1=TF[,sameSample]
immuneGene1=immuneGene[,sameSample]

#?????Լ???
outTab=data.frame()
for(i in row.names(TF1)){
  if(sd(TF1[i,])>1){
    for(j in row.names(immuneGene1)){
      x=as.numeric(TF1[i,])
      y=as.numeric(immuneGene1[j,])
      corT=cor.test(x,y)
      cor=corT$estimate
      pvalue=corT$p.value
      if((cor>corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(TF=i,immuneGene=j,cor,pvalue,Regulation="postive"))
      }
      if((cor< -corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(TF=i,immuneGene=j,cor,pvalue,Regulation="negative"))
      }
    }
  }
}
write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F)    

immuneGeneSig = read.table("tcgaUniCox.txt", row.names=1 ,header=T,sep="\t",check.names=F)
rownames(immuneGeneSig)=gsub("\\|","\\-",rownames(immuneGeneSig))
immuneGeneUp=immuneGeneSig[immuneGeneSig$HR>1,]
immuneGeneDown=immuneGeneSig[immuneGeneSig$HR<1,]
TFLabel=cbind(rownames(TF),"TF")
immuneGeneupLabel=cbind(rownames(immuneGeneUp),"highRiskIRG")
immuneGenedownLabel=cbind(rownames(immuneGeneDown),"lowRiskIRG")
nodeLabel=rbind(c("ID","type"),TFLabel,immuneGeneupLabel,immuneGenedownLabel)
write.table(nodeLabel,file="nodeType.txt",sep="\t",quote=F,col.names=F,row.names=F)