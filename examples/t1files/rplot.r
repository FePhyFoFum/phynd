args <- commandArgs(trailingOnly = TRUE)
a = read.table(args[1],sep=" ",h=T)
library(gplots)
png(args[2])
heatmap.2(as.matrix(a),dendrogram = 'none',sepwidth=c(0.05,0.05), sepcolor="white",
          lhei = c(0.2,1),margins=c(9,10),main=args[3],key=F,colsep=c(1:10000),rowsep=(1:10000), 
          Rowv=FALSE,Colv=FALSE,trace='none', cexRow=1.0,cexCol = 1,labRow="",
          breaks=seq(-1,1,length=4),col=c("red","white","blue"))
dev.off()
