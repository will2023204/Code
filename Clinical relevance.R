
library(plyr)
library(ggplot2)
library(ggpubr)
scoreFile="tcgaRisk.txt"   
cliFile="clinical.txt"         
trait="M"                   
setwd("..") 
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(score), row.names(cli))
rt=cbind(score[sameSample,], cli[sameSample,])
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
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
pdf(file=paste0(trait,".pdf"), width=3, height=4)
print(p)
dev.off()

rt2=rt[,c(trait, "riskScore")]
colnames(rt2)=c("trait", "riskScore")
type=levels(factor(rt2[,"trait"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
boxplot=ggboxplot(rt2, x="trait", y="riskScore", fill="trait",
		          xlab=trait,
		          ylab="riskScore",
		          legend.title=trait,
		          palette=bioCol
		          )+ 
	    stat_compare_means(comparisons=my_comparisons)
pdf(file="boxplot.pdf",width=4,height=4.5)
print(boxplot)
dev.off()

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
riskFile="tcgarisk.txt"      #?????ļ?
cliFile="clinical.txt"        #?ٴ??????ļ?
setwd("..")

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[order(risk$riskScore),] 
cli=read.table(cliFile,sep="\t",header=T,check.names=F,row.names=1)
samSample=intersect(row.names(risk), row.names(cli))
risk=risk[samSample,c("risk","riskScore"),drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(cli, risk)
write.table(rt,"clinical-risk.txt",quote = F,sep = "\t",row.names = T,col.names = T)

sigVec=c()
for(clinical in colnames(rt[,1:(ncol(rt)-1)])){
  data=rt[c("risk", clinical)]
  colnames(data)=c("risk", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  tableStat=table(data)
  stat=chisq.test(tableStat)
  pvalue=stat$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  sigVec=c(sigVec, paste0(clinical, Sig))
  print(tableStat)
  print(paste(clinical, pvalue, Sig, sep="\t"))
}
sigVec=c(sigVec, "risk")
colnames(rt)=sigVec
rt$risk=factor(rt$risk, levels=c("low","high"))

#rt=rt[apply(rt,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
bioCol=c("#0066FF","#FF9900","#FF0000","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
         "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
colorList=list()
j=0
for(cli in colnames(rt[,1:(ncol(rt)-1)])){
  cliLength=length(levels(factor(rt[,cli])))
  cliCol=bioCol[(j+1):(j+cliLength)]
  j=j+cliLength
  names(cliCol)=levels(factor(rt[,cli]))
  cliCol["unknow"]="grey75"
  colorList[[cli]]=cliCol
}
colorList[["risk"]]=c("low"="blue", "high"="red")

ha=HeatmapAnnotation(df=rt, col=colorList)
zero_row_mat=matrix(nrow=0, ncol=nrow(rt))
Hm=Heatmap(zero_row_mat, top_annotation=ha)

pdf(file="heatmap.pdf", width=7, height=5)
draw(Hm, merge_legend=TRUE, heatmap_legend_side="bottom", annotation_legend_side="bottom")
dev.off()



options(stringsAsFactors=F)
library(limma)
library(ggpubr)


expFile="expression.txt"         #?????????ļ?
cliFile="clinical.txt"          
setwd("..")     

rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#gene=read.table(geneFile,header=F,sep="\t",check.names=F)
#data=data[as.vector(gene[,1]),]

#group=sapply(strsplit(colnames(data),"\\-"),"[",4)
#group=sapply(strsplit(group,""),"[",1)
#group=gsub("2","1",group)
#data=data[,group==0]
#colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(data))
data=t(data)
data=avereps(data)

cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1) 
xlabel=colnames(cli)[1]

sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
rt=cbind(as.data.frame(data),cli)
write.table(rt,"clinical-gene expression.txt",quote = F,sep = "\t",row.names = T,col.names = T)

clinical="cli"
for(gene in colnames(rt)[1:(ncol(rt)-1)]){
  data=rt[c(gene,clinical)]
  colnames(data)=c("gene","clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  
  group=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  boxplot=ggboxplot(data, x="clinical", y="gene", color="clinical",
                    xlab=xlabel,
                    ylab=paste(gene,"expression"),
                    legend.title=xlabel,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  pdf(file=paste0(gene,".pdf"),width=5.5,height=5)
  print(boxplot)
  dev.off()
}
