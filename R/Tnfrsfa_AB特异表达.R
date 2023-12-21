library(ggplot2)
setwd("C:\\Users\\lcc\\Desktop\\ear")
data <- read.csv("Tnfrsfa_AB特异表达.csv",header = TRUE)

data$Cell_type <- factor(data$Cell_type,
                            levels=c('HC',
                                     'SC','mantle cells, anterior-posterior (A/P) cells and polar cells'), ordered=TRUE)    ####红色区域进行输入排列


ggplot(data,aes(x=name,y=logFC,fill=factor(Cell_type)))+ 
  geom_bar(position="dodge",stat="identity")+labs(fill="Cell_type")
  labs(x="gene",y="logFC")+theme(axis.text.x = element_blank(),axis.title.y = element_blank(),
                                axis.ticks.x = element_blank(),axis.title.x=element_blank())



