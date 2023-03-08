library(limma)
library(sva)
setwd("..")

rt = read.table("symbol.txt",header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
DDrepair=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
DDrepair=avereps(DDrepair)
genes = read.table("genes.txt",header=T,sep="\t",check.names=F)
abc = intersect(rownames(DDrepair),genes[,1])
DDrepair1=DDrepair[abc,]
write.table(DDrepair1,file="tcgaDNArepairExp.txt",sep="\t",quote=F,col.names=T)

rt = read.table("GEOmatrix.txt",header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo=avereps(geo)

sameGene=intersect(row.names(DDrepair1),row.names(geo))
DDrepairOut=DDrepair1[sameGene,]
geoOut=geo[sameGene,]

all=cbind(DDrepairOut,geoOut)
batchType=c(rep(1,ncol(DDrepairOut)),rep(2,ncol(geoOut)))
outTab=ComBat(all, batchType,par.prior=TRUE)
DDrepairOut=outTab[,colnames(DDrepairOut)]
DDrepairOut=rbind(ID=colnames(DDrepairOut),DDrepairOut)
write.table(DDrepairOut,file="tcga-DDR-Exp.share.txt",sep="\t",quote=F,col.names=F)
geoOut=outTab[,colnames(geoOut)]
geoOut=rbind(ID=colnames(geoOut),geoOut)
write.table(geoOut,file="geo-DDR-Exp.share.txt",sep="\t",quote=F,col.names=F)
