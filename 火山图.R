#### 火山图
library(ggpubr)
library(ggthemes)
res$logP <- -log10(res$padj)
res2<- data.frame(res)
### 分组
res$Group ="no significan"
res$Group[which((res$padj<0.05)& (res$log2FoldChange > 2))] ="up regulated"
res$Group[which((res$padj<0.05)& (res$log2FoldChange < -2))] ="down regulated"
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
res2 <-data.frame(res)
ggscatter(res2,x="log2FoldChange",y="logP",color = "Group",
          palette = c("#2f5688","#BBBBBB","#CC0000"),
          size = 1, label = res$Lable,font.label = 8,
          repel = T,xlab = "logFC",ylab = "-log(p_adj)")+theme_base()+
  geom_hline(yintercept = 1.3,linetype="dashed")+
  geom_vline(xintercept = c(-2,2),linetype="dashed")
