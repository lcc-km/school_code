library(edgeR)
#library(Glimma)
#BiocManager::install("Glimma")

#require("RColorBrewer"); display.brewer.all()
setwd("C:\\Users\\lucc\\Desktop\\微创文件\\9.ERA\\edgeR")
spe <- "ERA_"   ## 设置比较的物种

data <- read.csv("C:\\Users\\lucc\\Desktop\\微创文件\\9.ERA\\gene_count.csv", header=T, row.names=1)
keep <- rowSums(cpm(data)>2)>=2      #more than 2 counts per million in at least 2 repeats group
d <- data[keep,]
d_f <- d[,grep("F.*LHplus[0-9]",colnames(d))] 
###设置分组

targets <- data.frame(seq=c(1:24),
           sample=colnames(d_f),
           control=c("LHplus2","LHplus2","LHplus7","LHplus2",
                     "LHplus7","LHplus2","LHplus7","LHplus2",
                     "LHplus7","LHplus2","LHplus7","LHplus2",
                     "LHplus7","LHplus2","LHplus7","LHplus2",
                     "LHplus7","LHplus7","LHplus2","LHplus7",
                     "LHplus2","LHplus7","LHplus2","LHplus7"))



############ 计算表达差异 gene ##########################
###  normalize by TMM and save data to exp_study

y = DGEList(counts=d_f, group=targets$control)           ###6个样完整时用

y = calcNormFactors(y,method = "TMM") #get TMM nomalization factor # 计算样本内标准化因子
y$samples
y
colnames(y) <- targets$sample ###6个样完整时用
dim(y)
###  PCA analysis , 评估样品整体表达模式之间的相似程度
pp<-plotMDS(y,method = 'BCV',col=as.numeric(factor(targets$control)))
###  estimating the dispersion  or  Estimating dispersion
y1 <- y #get data
y2= estimateCommonDisp(y1,verbose=TRUE)    #计算普通的离散度
y2 = estimateTagwiseDisp(y2)   #计算基因间范围内的离散度
y2
plotBCV(y2)

##set threshold  of DE gene #############################################################
thrd <- 0.05  #set threshold of difference expression
chg_rd=1 ##set threshold of log2 change fold

#### Testing for DE genes
et <- exactTest(y2)     # 进行精确检验
top<-topTags(et)   # 输出排名靠前的差异表达基因信息
top
et.d <-data.frame(topTags(et,n=nrow(et$table)))
data0 <- et.d[abs(et.d[,"logFC"])>chg_rd & et.d[,"FDR"]<thrd,]
filename <- paste(spe,"true","exact_FC.txt",sep ="")
write.table(data0,filename,quote=F,sep = "\t")

summary(de <- decideTestsDGE(et)); # 默认选取FDR = 0.05为阈值
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")

ggplot(res,aes(x= log2FoldChange,
               y= -1*log10(res$pval),colour=threshold))+
  xlab("log2 fold-change")+ylab("-log10 p-value")+geom_point() 



