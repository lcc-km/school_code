library(edgeR)
require("RColorBrewer"); display.brewer.all()
setwd("C:\\Users\\lcc\\Desktop\\ear")
spe <- "ear_"   ## 设置比较的物种

data <-  read.csv("gene_count.csv", header=T, row.names=1)
keep <- rowSums(cpm(data)>2)>=2      #more than 2 counts per million in at least 2 repeats group
d <- data[keep,]

###设置分组

targets <- data.frame(seq=c(1:6),
                      control=c("AB","AB","AB","Tnf","Tnf","Tnf"),
                      Label=c("AB-1","AB-1-2","AB-2","Tnf-1","Tnf-1-2","Tnf-2"))    ###6个样完整时用

############ 计算表达差异 gene ##########################
###  normalize by TMM and save data to exp_study

y = DGEList(counts=d, group=targets$control)           ###6个样完整时用

y = calcNormFactors(y,method = "TMM") #get TMM nomalization factor # 计算样本内标准化因子
y$samples
y
colnames(y) <- targets$Label   ###6个样完整时用
dim(y)
###  PCA analysis , 评估样品整体表达模式之间的相似程度


pp<-plotMDS(y,method="bcv",col=c(rep("royalblue",3),rep("red3",3)))
#################t-sne###########
plot(tsne_out$Y,col=c(rep("royalblue",3),rep("red3",3)))
###  estimating the dispersion  or  Estimating dispersion
y1 <- y #get data
y2= estimateCommonDisp(y1,verbose=TRUE)    #计算普通的离散度
y2 = estimateTagwiseDisp(y2)   #计算基因间范围内的离散度
y2
plotBCV(y2)

#### Testing for DE genes
et <- exactTest(y2)     # 进行精确检验
top<-topTags(et)   # 输出排名靠前的差异表达基因信息
top
et.d <-data.frame(topTags(et,n=nrow(et$table)))
filename <- paste(spe,"Tnf-AB","exact_FC.txt",sep ="")
write.table(et.d,filename,quote=F,sep = "\t")

summary(de <- decideTestsDGE(et)); # 默认选取FDR = 0.05为阈值
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")

