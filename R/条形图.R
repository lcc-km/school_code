library(ggplot2)
setwd("C:\\Users\\lcc\\Desktop\\inject")
data<-read.csv("数据统计.csv",header=TRUE)
data$percent<-paste((substr(data$count/data$sun*100,1,4)),"%",sep="")
data$process<-factor(data$process,levels=c("mock","Enhancer"),ordered=TRUE)
data_12<-data[data$period=="12h",]
data_24<-data[data$period=="24h",]
data_36<-data[data$period=="36h",]
data_48<-data[data$period=="48h",]
ggplot(data_12,aes(x=process,y=count,fill=factor(light)))+
  geom_bar(position="dodge",stat="identity")+labs(fill="light")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="top")+
  geom_text(aes(label =percent ), size = 2.5,position = position_dodge(0.9))+
  geom_text(aes(label =light,y=3), size = 3,position = position_dodge(0.9))+
  guides(fill=FALSE)

ggplot(data_24,aes(x=process,y=count,fill=factor(light)))+
  geom_bar(position="dodge",stat="identity")+labs(fill="light")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="top")+
  geom_text(aes(label =percent ), size = 2.5,position = position_dodge(0.9))+
  geom_text(aes(label =light,y=3), size = 3,position = position_dodge(0.9))+
  guides(fill=FALSE)

ggplot(data_36,aes(x=process,y=count,fill=factor(light)))+
  geom_bar(position="dodge",stat="identity")+labs(fill="light")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="top")+
  geom_text(aes(label =percent ), size = 2.5,position = position_dodge(0.9))+
  geom_text(aes(label =light,y=3), size = 3,position = position_dodge(0.9))+
  guides(fill=FALSE)

ggplot(data_48,aes(x=process,y=count,fill=factor(light)))+
  geom_bar(position="dodge",stat="identity")+labs(fill="light")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="top")+
  geom_text(aes(label =percent ), size = 2.5,position = position_dodge(0.9))+
  geom_text(aes(label =light,y=3), size = 3,position = position_dodge(0.9))+
  guides(fill=FALSE)

