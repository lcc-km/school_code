setwd("G:\\backup\\上班\\转录组数据\\自己处理\\conA_wt_VS-conA-ko")
library(dplyr)
library(ggplot2)
d1<-read.table("mm_conA-KO-SP-VS-conA-WT-SP_exact_FC1_FDR0.05.txt")
d2<-read.table("mm_conA-KO-BM-VS-conA-WT-BM_exact_FC1_FDR0.05.txt")
d1$gene <- row.names(d1)
d2$gene <- row.names(d2)
d1$lable <-"SP"
d2$lable <-"BM"
d1_d2<-rbind(d1,d2)

d_low <-d1_d2[d1_d2$logFC<0,]
d_high <-d1_d2[d1_d2$logFC>0,]


#ko_up_gene<- d[d$gene%in%low,]
#ko_down_gene<- d[d$gene%in%high,]

#write.table(low,"ko_up_gene",quote=F,sep = "\t",row.names = F,col.names = "gene")
#write.table(high,"ko_down_gene",quote=F,sep = "\t",row.names = F,col.names = "gene")

write.table(d_high,"conA-KO-VS-conA-WT_up_gene_data.txt",quote=F,sep = "\t",row.names = F)
write.table(d_low,"conA-KO-VS-conA-WT_down_gene_data.txt",quote=F,sep = "\t",row.names = F)

#####go
data_BP<-  read.csv("down_go.txt", header=T,sep = "\t")
colnames(data_BP)<-c("GO_biological_process","ref","count","expected",
                     "+-","fold_Enrichment","P_value","FDR")
data_BP$GO_biological_process <- sapply(strsplit(as.character(data_BP$GO_biological_process),"(",fixed= T), "[", 1)

data_BP2<-arrange(data_BP,count,FDR)

ggplot(data_BP2,aes(fold_Enrichment,GO_biological_process))+ 
  geom_point(aes(size=count))+
  geom_point(aes(size=count,color=-1*log10(FDR)))+
  scale_color_gradient(low="#0072B2",high = "#CC79A7")+
  labs(x="fold Enrichment",y="GO biological process",
       size="Gene number",color="-log10(FDR)")

ggplot(data_BP2,aes(-1*log10(FDR),GO_biological_process))+
  geom_bar(stat="identity",fill = "#CC79A7")+ 
  labs(x="-log10(FDR)",y="GO biological process")+
  theme_classic()+geom_text(aes(label=count),hjust = -1)


####kegg
#down
data_kegg_down<-  read.csv("down-kegg.txt", header=T,sep = "\t")
data_kegg_down<-data_kegg_down%>%arrange(Corrected.P.Value,P.Value)
data_kegg_down<-data_kegg_down[data_kegg_down$Corrected.P.Value<0.05,]

ggplot(data_kegg_down,aes(-1*log10(Corrected.P.Value),Term))+
  geom_bar(stat="identity",fill = "#CC79A7",width=0.3)+
  labs(x="-log10(FDR)",y="kegg term")+
  theme_classic()+geom_text(aes(label=Input.number),hjust = -1)

#up
data_kegg_up<-  read.csv("up-kegg.txt", header=T,sep = "\t")
data_kegg_up<-data_kegg_up%>%arrange(Corrected.P.Value,P.Value)
data_kegg_up<-data_kegg_up[data_kegg_down$Corrected.P.Value<0.05,]

ggplot(data_kegg_up,aes(-1*log10(Corrected.P.Value),Term))+
  geom_bar(stat="identity",fill = "#CC79A7",width=0.3)+
  labs(x="-log10(FDR)",y="kegg term")+
  theme_classic()+geom_text(aes(label=Input.number),hjust = -1)



