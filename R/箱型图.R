library(ggplot2)
setwd("D:\\backup\\our_data\\RNA\\VQSR\\10\\gene_area\\expression_quantity")
data<- read.table("RNA-seq_data_count_2.txt",header=TRUE,sep="\t")
#data$type=factor(data$type,levels=c("total SNP","aftet VQSR","aftet VQSR number of heterozygosis","aftet VQSR number of heterozygosis locate in gene"))
ggplot(data,aes(x=class,y=number))+
  geom_boxplot(aes(fill=type))+theme(legend.position="top")

