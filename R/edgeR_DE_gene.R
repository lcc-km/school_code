library(edgeR)
setwd("C:\\Users\\lcc\\Desktop\\our_data\\RNA")
spe <- "zeb_"   ## 设置比较的物种

data <-  read.csv("gene_count.csv", header=T, row.names=1)
keep <- rowSums(cpm(data)>2)>=2      #more than 2 counts per million in at least 2 repeats group
d <- data[keep,]

###设置分组
#group <- factor(gsub("\..","",colnames(d)))
#targets <- data.frame(seq=c(1:36),
#                      control=c(c(rep("A",6),rep("B",6),rep("C",6),
#                                  rep("D",6),rep("E",6),rep("F",6))),
#                      Label=c(paste("A", 1:6, sep=""),
#                              paste("B", 1:6, sep=""),
#                              paste("C", 1:6, sep=""),
#                              paste("D", 1:6, sep=""),
#                              paste("E", 1:6, sep=""),
#                              paste("F", 1:6, sep="")))    ###36个样完整时用

group <- c(1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5,6,6,6,6,6,6)
#substring(colnames(d),first=7)
##conditions <- factor(colnames(d))
##conditions <- factor(substring(colnames(d),first=7,last=nchar(colnames(d))-1)) #get col.names and factor as conditions

############ 计算表达差异 gene ##########################
###  normalize by TMM and save data to exp_study
y = DGEList(counts=d, group=group)
###### y = DGEList(counts=d, group=targets$control)           ###36个样完整时用


y = calcNormFactors(y,method = "TMM") #get TMM nomalization factor
y$samples
y
###colnames(y) <- targets$Label    ###36个样完整时用
dim(y)
###  PCA analysis , 评估样品整体表达模式之间的相似程度
plotMDS(y,method="bcv",col="firebrick")

### 样品没有重复，人为设置BVC
#BVC=0.1


###  estimating the dispersion  or  Estimating dispersion
y1 <- y #get data
y2= estimateCommonDisp(y1,verbose=TRUE)
#y2$common.dispersion=BVC
y2 = estimateTagwiseDisp(y2)
y2
plotBCV(y2)

#### Testing for DE genes

et <- exactTest(y2)

write.table(et$table,"all_in_one.txt",quote=F,sep = "\t")
###
t1 <- c("A","B","C","D","E","F")
t1_combn <-t(combn(t1,2))


##set threshold  of DE gene #############################################################
thrd <- 0.05  #set threshold of difference expression
chg_rd=1 ##set threshold of log2 change fold
##set threshold#########################################################################

## 比较差异gene

for(i in 1:length(t)){
  test <- t[[i]] ##get names of test pairs
  comp <- paste(test[2],test[1],sep = "-")   ## test[2]-test[1]
  comp

###   exactTest & write file
  et0 <- exactTest(y2,pair = test)   ###
  topTags(et0)
  et0.d <-data.frame(topTags(et0,n=nrow(et0$table)))
  filename1 <- paste(spe,comp,"_all_exact_FC.txt",sep ="")
  write.table(et0.d,filename1,quote=F,sep = "\t")


### drawing DE gene plot (optioinal)
  de <- decideTestsDGE(et0,p=thrd,lfc = chg_rd)
  summary(de) ##get number of different exp gene, -1 down regulate#,1 up regulate, 0 no difference
  plotSmear(et0,de.tags = rownames(y2)[as.logical(de)], pch=20, cex=0.6)
  abline(h=c(0-chg_rd,chg_rd),col="blue")

  filename <- paste(spe,comp,"_exact_FC",chg_rd,"_FDR",thrd,".pdf",sep ="")
  pdf(file = filename,width = 5,height = 5,title = comp)
  plotSmear(et0,de.tags = rownames(y2)[as.logical(de)])
  abline(abline(h = c(-1, 0, 1), col = c("dodgerblue", "yellow", "dodgerblue")))
  title(paste(spe,comp))
  dev.off()
}
