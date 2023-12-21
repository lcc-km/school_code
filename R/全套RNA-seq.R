library(edgeR)
library(ggpubr)
library(ggthemes)
library(tidyr)
library(pheatmap)
library(stringr)
require("RColorBrewer"); display.brewer.all()
setwd("D:\\学校\\ear")
spe <- "mm_"   ## 设置比较的物种
gene_count_file <- paste('gene_count.csv',sep = "\\")
data <- read.csv(gene_count_file, header=T, row.names=1)
keep <- rowSums(cpm(data)>2)>=2      #more than 2 counts per million in at least 2 repeats group
d <- data[keep,]

###设置分组

targets <- data.frame(seq=c(1:6),
                      control=c("AB","AB","AB","Tnf","Tnf","Tnf"),
                      Label=c("AB-1","AB-1-2","AB-2","Tnf-1","Tnf-1-2","Tnf-2"))  

############ 计算表达差异 gene ##########################
###  normalize by TMM and save data to exp_study

y = DGEList(counts=d, group=targets$control)         

y = calcNormFactors(y,method = "TMM") #get TMM nomalization factor
y$samples
y

colnames(y) <- targets$sample
dim(y)
###  PCA analysis , 评估样品整体表达模式之间的相似程度
#glMDPlot(y,labels=targets2$control)
pp<-plotMDS(y,method = 'BCV',col=as.numeric(factor(targets$control)))
pp<-plotMDS(y,col=as.numeric(factor(targets$control)))
###  estimating the dispersion  or  Estimating dispersion
y1 <- y #get data
y2= estimateCommonDisp(y1,verbose=TRUE)
y2 = estimateTagwiseDisp(y2)
y2
plotBCV(y2)

#### Testing for DE genes
et <- exactTest(y2)

##设置比较组
t1 <- unique(factor(targets$control))
t1_combn <-t(combn(t1,2))

t1_combn
##set threshold  of DE gene #############################################################
thrd <- 0.01  #set threshold of difference expression
chg_rd=1 ##set threshold of log2 change fold
##set threshold#########################################################################

## 比较差异gene

for(i in 1:length(t1_combn[,1])) {
  test <- t1_combn[i,] ##get names of test pairs
  comp <- paste(test[2],test[1],sep = "-")   ## test[2]-test[1]
  comp
  
  ###   exactTest & write file
  et0 <- exactTest(y2,pair = test)
  #print(summary(de <- decideTestsDGE(et0))) 
  topTags(et0)
  et0.d <-data.frame(topTags(et0,n=nrow(et0$table)))
  filename1 <- paste(spe,comp,"_all_exact_FC.txt",sep ="")
  write.table(et0.d,filename1,quote=F,sep = "\t")
  
  
  # drawing DE gene plot (optioinal)
  de <- decideTestsDGE(et0,p=thrd,lfc = chg_rd)
  
  print(summary(de) ) ##get number of different exp gene, -1 down regulate??1 up regulate, 0 no difference
  de_inf <- summary(de) 
  filename_de<-paste(spe,comp,"_differential_gene_number.txt",sep ="")
  write.table(de_inf,filename_de,quote=F,sep = "\t")
  ### 火山图
  res <- et0.d 
  colnames(res) <- c("log2FoldChange","logCPM","pvalue","padj")
  res$logP <- -log10(res$padj)
  ### 分组
  res$Group ="no significan"
  res$Group[which((res$padj<thrd)& (res$log2FoldChange > chg_rd))] ="up regulated"
  res$Group[which((res$padj<thrd)& (res$log2FoldChange < -chg_rd))] ="down regulated"
  table(res$Group)
  ### 提取 top10 gene
  res$Lable = ""
  res <- res[order(res$padj),]
  res$ID <- rownames(res)
  up_genes <- head(res$ID[which(res$Group == "up regulated")],10)
  down_genes <- head(res$ID[which(res$Group == "down regulated")],10)
  top_10_gene <- c(as.character(up_genes),as.character(down_genes))
  res$Lable[match(top_10_gene,res$ID)] <- top_10_gene
  ### plot
  p<-ggscatter(res,x="log2FoldChange",y="logP",color = "Group",
               palette = c("#2f5688","#BBBBBB","#CC0000"),
               size = 1, 
               #label = res$Lable,font.label = 8,
               repel = T,xlab = "logFC",ylab = "-log(p_adj)")+theme_base()+
    labs(title=paste(spe,comp))+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_hline(yintercept = -log10(thrd) ,linetype="dashed")+
    geom_vline(xintercept = c(-2^chg_rd,2^chg_rd),linetype="dashed")
  p
  filename <- paste(spe,comp,"_exact_FC",chg_rd,"_FDR",thrd,"_volcano",".pdf",sep ="")
  ggsave(filename,p,width=8, height=6) # 可以指定大小，如宽为12cm，高为10cm，需要指定保存路径
  
  
  #get DE gene text file
  data0 <- et0.d[abs(et0.d[,"logFC"])>chg_rd & et0.d[,"FDR"]<thrd,]
  filename <- paste(spe,comp,"_exact_FC",chg_rd,"_FDR",thrd,".txt",sep ="")
  
  write.table(data0,filename,quote=F,sep = "\t")
  
  ###### heatmap
  data1 <- data0
  data1 <- head(data1,n=20)
  data1$gene <- row.names(data1)
  d$gene <- row.names(d)
  h_d <- d[d$gene%in%data1$gene,]
  h_d2 <- h_d[,c(1:(length(colnames(h_d))-1))]
  h_d3 <- scale(h_d2)

  heatmap_name <- paste(spe,comp,"_exact_FC",chg_rd,"_FDR",thrd,"_heatmap.pdf",sep ="")
  #pdf(file = heatmap_name,width = 6,height = 6)
  p<-pheatmap(h_d3,border=FALSE)
  ggsave(heatmap_name,p,width=8, height=6) 
  #dev.off()
}






