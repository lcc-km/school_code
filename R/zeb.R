library("edgeR")
setwd("C:\\Users\\49248\\Desktop\\rna-seq\\zeb")
spe <- "zeb"   ## 设置比较的物种
foldChange=1
padj=0.05
options(stringsAsFactors = FALSE)   ##遇到字符串之后，不将其转换为factors，仍然保留为字符串格式
data <-  read.csv("gene_count.csv", header=TRUE, row.names=1)
data <- as.matrix(data)       ##转化为矩阵
dim(data)   ##检查一下数量
keep <- rowSums(cpm(data)>2)>=2     ##在至少2个重复组中每百万次数超过2次
d <- data[keep,] 
d <- as.matrix(d)
dim(d)     ##检查筛选过后的数量
group <- factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8))  #按样本类型分类
y <- DGEList(counts=d,group=group)
colnames(y)
y = calcNormFactors(y,method = "TMM")
y
plotMDS(y)
y <- estimateCommonDisp(y, verbose=TRUE)
plotSmear(y)
y <- estimateTagwiseDisp(y)
plotBCV(y)

et <- exactTest(y)
top <- topTags(et)
top
summary(de <- decideTestsDGE(et))
detags <- rownames(y)[as.logical(de)];
plotSmear(et, de.tags=detags);
abline(h=c(-1, 1), col="blue");

ordered_tags <- topTags(et, n=100000)
allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts
write.table(diff,file="edgerOut.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)
normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)
heatmapData <- newData[rownames(diffSig),]

#volcano
pdf(file="vol.pdf") #注意输出的火山图的路径#
xMax=max(-log10(allDiff$FDR))+1
yMax=12
plot(-log10(allDiff$FDR), allDiff$logFC, xlab="-log10(FDR)",ylab="logFC",
     main="Volcano",yaxs="i",pch=20, cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC>foldChange,]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC<(-foldChange),]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()

#heatmap

hmExp=log10(heatmapData+0.001)
library('gplots')
hmMat=as.matrix(hmExp)
pdf(file="heatmap.pdf",width=60,height=90)#注意输出的热图的储存路径#
par(oma=c(10,3,3,7))
heatmap.2(hmMat,col='greenred',trace="none")
dev.off()

b <- sort(data$logFC,decreasing= T)


setwd("C:\\Users\\49248\\Desktop\\parallel comparisonRNA-seq\\28-zeb")
data <-  read.csv("zebGills-Brain_exact_FC1_FDR0.05.txt", header=TRUE ,sep="\t")