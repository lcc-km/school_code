library(ggplot2)
setwd("C:\\Users\\pro3\\Desktop\\ear")
data <- read.csv("gene.GO_edgeR.txt",header=TRUE,sep = "\t")
ID <- read.csv("TNF_AB_edgeR_lucc.txt",header=TRUE,sep = "\t")
data_work <- read.csv("data_work_1.txt",header=TRUE,sep = "\t")
data1<-merge(data,data_work,all=FALSE)
data2 <-merge(data1,ID,by.x = "name",by.y = "gene_name")
write.table(data2,"go_enrichment_edgeR.csv",sep = ",")
data2$Description <- factor(data2$Description, 
                            levels=c('positive regulation of Notch signaling pathway',
                                                        'regulation of Notch signaling pathway',
                                                        'Notch signaling pathway',
                                                        'negative regulation of Notch signaling pathway',
                                                        'Notch binding',
                                                        'non-canonical Wnt signaling pathway',
                                                        'regulation of non-canonical Wnt signaling pathway',
                                                        'negative regulation of canonical Wnt signaling pathway',
                                                        'negative regulation of Wnt signaling pathway',
                                                        'canonical Wnt signaling pathway',
                                                        'regulation of canonical Wnt signaling pathway',
                                                        'regulation of Wnt signaling pathway',
                                                        'Wnt signaling pathway',
                                                        'cell-cell signaling by wnt',
                                                        'positive regulation of Wnt signaling pathway',
                                                        'positive regulation of canonical Wnt signaling pathway',
                                                        'Wnt-protein binding',
                                                        'Wnt-activated receptor activity',
                                                        'I-kappaB kinase/NF-kappaB signaling',
                                                        'positive regulation of I-kappaB kinase/NF-kappaB signaling',
                                                        'NIK/NF-kappaB signaling',
                                                        'regulation of I-kappaB kinase/NF-kappaB signaling',
                                                        'regulation of NIK/NF-kappaB signaling',
                                                        'positive regulation of JNK cascade',
                                                        'regulation of JNK cascade',
                                                        'JNK cascade'), ordered=TRUE)

ggplot(data2, aes(x=name, y=Description, fill=logFC)) + 
  geom_tile()+
  scale_fill_gradient2(name="logFC", low = "cyan",high = "orange",mid = "white")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))



