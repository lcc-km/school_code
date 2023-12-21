library(edgeR)
require("RColorBrewer"); display.brewer.all()
setwd("C:\\Users\\lcc\\Desktop\\our_data\\RNA\\raw")
spe <- "zeb_"   ## 设置比较的物种

data <-  read.csv("gene_count.raw.csv", header=T, row.names=1)
keep <- rowSums(cpm(data)>2)>=2      #more than 2 counts per million in at least 2 repeats group
d <- data[keep,]

###设置分组

targets <- data.frame(seq=c(1:36),
                      control=c(c(rep("Af",3),rep("Am",3),rep("Bf",3),
                                  rep("Bm",3),rep("Cf",3),rep("Cm",3),
                                  rep("Df",3),rep("Dm",3),rep("Ef",3),
                                  rep("Em",3),rep("Ff",3),rep("Fm",3))),
                      Label=c(paste("Af", 1:3, sep=""),
                              paste("Am", 1:3, sep=""),
                              paste("Bf", 1:3, sep=""),
                              paste("Bm", 1:3, sep=""),
                              paste("Cf", 1:3, sep=""),
                              paste("Cm", 1:3, sep=""),
                              paste("Df", 1:3, sep=""),
                              paste("Dm", 1:3, sep=""),
                              paste("Ef", 1:3, sep=""),
                              paste("Em", 1:3, sep=""),
                              paste("Ff", 1:3, sep=""),
                              paste("Fm", 1:3, sep="")))    ###36个样完整时用

############ 计算表达差异 gene ##########################
###  normalize by TMM and save data to exp_study

y = DGEList(counts=d, group=targets$control)           ###36个样完整时用

y = calcNormFactors(y,method = "TMM") #get TMM nomalization factor
y$samples
y
colnames(y) <- targets$Label   ###36个样完整时用
dim(y)
###  PCA analysis , 评估样品整体表达模式之间的相似程度


pp<-plotMDS(y,method="bcv",col=c(rep("royalblue",3),rep("red3",3),rep("forestgreen",3),
                                 rep("gold",3),rep("black",3),rep("pink",3),
                                 rep("blue",3),rep("brown1",3),rep("coral2",3),
                                 rep("azure4",3),rep("aquamarine4",3),
                                 rep("aquamarine4",3)))
#################t-sne###########
library(Rtsne)
tsne_data <- t(y$counts)
colours <- factor(tsne_data,levels=c("Af","Am","Bf",
                                     "Bm","Cf","Cm",
                                     "Df","Dm","Ef",
                                     "Em","Ff","Fm"))
colours <- colorRampPalette(c("royalblue", "red3",
                              "forestgreen", "gold",
                              "black","pink","blue","brown1",
                              "coral2",
                              "azure4","aquamarine4","darkgray"))(length(unique(colours)))[factor(colours)]
tsne_out <- Rtsne(
  tsne_data,
  dims = 2,
  pca = FALSE,
  perplexity = 11,
  theta = 0.5,
  max_iter = 5000,
  num_threads = 8
)

plot(tsne_out$Y,col=c(rep("royalblue",3),rep("red3",3),rep("forestgreen",3),
                      rep("gold",3),rep("black",3),rep("pink",3),
                      rep("blue",3),rep("brown1",3),rep("coral2",3),
                      rep("azure4",3),rep("aquamarine4",3),
                      rep("orange2",3)))
###  estimating the dispersion  or  Estimating dispersion
y1 <- y #get data
y2= estimateCommonDisp(y1,verbose=TRUE)
y2 = estimateTagwiseDisp(y2)
y2
plotBCV(y2)

#### Testing for DE genes
et <- exactTest(y2)

##设置比较组
t1 <- c("A","B","C","D","E","F")
t1_combn <-t(combn(t1,2))

##set threshold  of DE gene #############################################################
thrd <- 0.05  #set threshold of difference expression
chg_rd=1 ##set threshold of log2 change fold
##set threshold#########################################################################

## 比较差异gene

for(i in 1:length(t1_combn[,1])) {
  test <- t1_combn[i,] ##get names of test pairs
  comp <- paste(test[2],test[1],sep = "-")   ## test[2]-test[1]
  comp

  ###   exactTest & write file
  et0 <- exactTest(exp_study1,pair = test,dispersion = BVC)
  topTags(et0)
  et0.d <-data.frame(topTags(et0,n=nrow(et0$table)))
  filename1 <- paste(spe,comp,"_all_exact_FC.txt",sep ="")
  write.table(et0.d,filename1,quote=F,sep = "\t")


  # drawing DE gene plot (optioinal)
  de <- decideTestsDGE(et0,p=thrd,lfc = chg_rd)
  summary(de) ##get number of different exp gene, -1 down regulate??1 up regulate, 0 no difference
  plotSmear(et0,de.tags = rownames(exp_study1)[as.logical(de)])
  abline(h=c(0-chg_rd,chg_rd),col="blue")

  filename <- paste(spe,comp,"_exact_FC",chg_rd,"_FDR",thrd,".pdf",sep ="")
  pdf(file = filename,width = 5,height = 5,title = comp)
  plotSmear(et0,de.tags = rownames(exp_study1)[as.logical(de)])
  abline(h=c(0-chg_rd,chg_rd),col="blue")
  title(paste(spe,comp))
  dev.off()


  #get DE gene text file
  data0 <- et0.d[abs(et0.d[,"logFC"])>chg_rd & et0.d[,"FDR"]<thrd,]
  filename <- paste(spe,comp,"_exact_FC",chg_rd,"_FDR",thrd,".txt",sep ="")
  write.table(data0,filename,quote=F,sep = "\t")

}
