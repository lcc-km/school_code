library(ggplot2)
library(dplyr)
library(tidyverse)
#install.packages("tidyverse")

Args<-commandArgs(T)
spe<-Args[1]
wd<-Args[2]


spe <- paste(spe,sep = "")
setwd(wd)
setwd("C:\\Users\\lucc\\Desktop\\微创文件\\9.ERA")
go <-read.table("47gene_GO.txt",sep = "\t",header = T)

colnames(go)


colnames(go) <- c("GO.biological.process","ref","count","expected","under","fold.Enrichment","raw.P.value","FDR")
go$GO_biological_process <- sapply(strsplit(as.character(go$GO.biological.process),"(",fixed= T), "[", 1)

go <- go[order(-go$count,go$FDR),]

top <- head(go,15)

p1 <-ggplot(top,aes(x=fold.Enrichment,y=reorder(GO_biological_process,fold.Enrichment)))+
  geom_point(aes(size=count,color=-1*log10(FDR)))+
  scale_color_gradient(low="#b0e0e6",high = "#Ff8c00")+
  labs(x="fold.Enrichment",y="GO biological process",title = "mature HC(cluster 16)",
       size="Gene number",color="-log10(FDR)")+theme_classic()  + 
  theme(axis.text.x = element_text(size=10),axis.text.y=element_text(size=10),
        axis.title.x= element_text(size=10),axis.title.y= element_text(size=10),
        plot.title = element_text(hjust = 0.5),title = element_text(size=10))#scale_fill_manual(values=c("#478bA2","#DDC1C6","#F2A490","#E9765B"))
p1
p1n <- paste(spe,"_GO_plot_1.pdf",sep = "")
ggsave(p1, file=p1n, width=9, height=4)

p2 <- ggplot(top,aes(x=-log10(FDR),y=reorder(GO_biological_process,-log10(FDR))))+
  geom_bar(stat="identity",fill = "#CC79A7")+ 
  labs(x="-log10(FDR)",y="GO biological process")+
  theme_classic()+geom_text(aes(label=count),hjust = -1)
p2
p2n <- paste(spe,"_GO_plot_2.pdf",sep = "")
ggsave(p2, file=p2n, width=10, height=4)


#########3 
kegg <-read.table("47gene_KEGG.txt",sep = "\t",header = T)
colnames(kegg) <- c( "Term","Database","ID","Input.number","Background.number",
                     "P.Value","FDR","Input","Hyperlink")
kegg <- kegg[order(-kegg$Input.number,kegg$FDR),]
top_kegg <- head(kegg,15)



p4 <- ggplot(top_kegg,aes(x=-log10(FDR),y=reorder(Term,-log10(FDR))))+
  geom_bar(stat="identity",fill = "#CC79A7")+ 
  labs(x="-log10(FDR)",y="GO biological process")+
  theme_classic()+geom_text(aes(label=Input.number),hjust = -1)
p4
p4n <- paste(spe,"_kegg_plot_1.pdf",sep = "")
ggsave(p4, file=p4n, width=10, height=4)
