  ##该函数再dplyr包存在的情况下会出错，选先删除运行了，再重新装???  detach("dplyr", unload=TRUE)
get_GO_data <- function(OrgDb, ont, keytype) {
  GO_Env <- get_GO_Env()
  use_cached <- FALSE
  
  if (exists("organism", envir=GO_Env, inherits=FALSE) &&
      exists("keytype", envir=GO_Env, inherits=FALSE)) {
    
    org <- get("organism", envir=GO_Env)
    kt <- get("keytype", envir=GO_Env)
    
    if (org == DOSE:::get_organism(OrgDb) &&
        keytype == kt &&
        exists("goAnno", envir=GO_Env, inherits=FALSE)) {
      ## https://github.com/GuangchuangYu/clusterProfiler/issues/182
      ## && exists("GO2TERM", envir=GO_Env, inherits=FALSE)){
      
      use_cached <- TRUE
    }
  }
  
  if (use_cached) {
    goAnno <- get("goAnno", envir=GO_Env)
  } else {
    OrgDb <- GOSemSim:::load_OrgDb(OrgDb)
    kt <- keytypes(OrgDb)
    if (! keytype %in% kt) {
      stop("keytype is not supported...")
    }
    
    kk <- keys(OrgDb, keytype=keytype)
    goAnno <- suppressMessages(
      select(OrgDb, keys=kk, keytype=keytype,
             columns=c("GOALL", "ONTOLOGYALL")))
    
    goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ])
    
    assign("goAnno", goAnno, envir=GO_Env)
    assign("keytype", keytype, envir=GO_Env)
    assign("organism", DOSE:::get_organism(OrgDb), envir=GO_Env)
  }
  
  if (ont == "ALL") {
    GO2GENE <- unique(goAnno[, c(2,1)])
  } else {
    GO2GENE <- unique(goAnno[goAnno$ONTOLOGYALL == ont, c(2,1)])
  }
  
  GO_DATA <- DOSE:::build_Anno(GO2GENE, get_GO2TERM_table())
  
  goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
  goOnt <- goOnt.df[,2]
  names(goOnt) <- goOnt.df[,1]
  assign("GO2ONT", goOnt, envir=GO_DATA)
  return(GO_DATA)
}
get_GO_Env <- function () {
  if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
  }
  get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
}
get_GO2TERM_table <- function() {
  GOTERM.df <- get_GOTERM()
  GOTERM.df[, c("go_id", "Term")] %>% unique
}
get_GOTERM <- function() {
  pos <- 1
  envir <- as.environment(pos)
  if (!exists(".GOTERM_Env", envir=envir)) {
    assign(".GOTERM_Env", new.env(), envir)
  }
  GOTERM_Env <- get(".GOTERM_Env", envir = envir)
  if (exists("GOTERM.df", envir = GOTERM_Env)) {
    GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
  } else {
    GOTERM.df <- toTable(GOTERM)
    assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
  }
  return(GOTERM.df)
}

GO_DATA <- get_GO_data("org.Dr.eg.db","ALL", "SYMBOL")

library(org.Dr.eg.db)

library(DOSE)
library(topGO)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(stringr)
library(GOplot)
library(enrichplot)
library(ggplot2)
library(cowplot)
library(dplyr) 
 
BiocManager::install("dplyr")

BiocManager::install("bitr")
##################### 
BiocManager::install(c("DOSE","topGO","clusterProfiler",
                   "pathview","stringr","GOplot",
                   "org.Dr.eg.db","enrichplot",
                     "cowplot","dplyr"))





setwd("D:\\backup\\enhancer paper")
#####raw data####
data_raw <-read.csv("183_gene.txt",header=F,sep="\t")
keep1 <- data_raw$FDR<=0.05
data <- data_raw[keep1,]
keep_high <- data$logFC>=0
keep_low <- data$logFC<=0
data_high <- data[keep_high,]
data_low <- data[keep_low,]
########
data_high <- read.csv("ear_RNA-seq_result_high_ncbi_ID.txt",header=FALSE,sep="\t")
data_high_FC <-data.frame(data_high$V1,data_high$V5)
data_low <- read.csv("ear_RNA-seq_result_low_ncbi_ID.txt",header=FALSE,sep="\t")
keytypes(org.Dr.eg.db)  #查询斑马鱼基因注释包中支持的ID类型
data <- rbind(data_high,data_low)
geneID <- bitr(data_raw$V1, fromType = "ENSEMBL",
       toType = c("ENSEMBL", "SYMBOL"),
       OrgDb = org.Dr.eg.db)
data2<-merge(data_raw,geneID,by.x="ID",by.y="ENSEMBL")

###GO enrichment
gene.GO_high_BP <- enrichGO(gene = data_raw$V1,
                            OrgDb = org.Dr.eg.db,
                            keyType = "ENSEMBL",
                            ont = "BP",
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            readable = T)
gene.GO_high_BP

gene.GO_high_MF <- enrichGO(gene = data_high$V1,
                            OrgDb = org.Dr.eg.db,
                            keyType = "ENSEMBL",
                            ont = "MF",
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            readable = T)
gene.GO_high_MF

gene.GO_high_CC <- enrichGO(gene = data_high$V1,
                            OrgDb = org.Dr.eg.db,
                            keyType = "ENSEMBL",
                            ont = "CC",
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            readable = T)
gene.GO_high_CC

#######all
gene.GO_BP <- enrichGO(gene = data$ID,
                            OrgDb = org.Dr.eg.db,
                            keyType = "ENSEMBL",
                            ont = "BP",
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 1,
                            qvalueCutoff = 1,
                            readable = T)
gene.GO_BP

gene.GO_MF <- enrichGO(gene = data$ID,
                            OrgDb = org.Dr.eg.db,
                            keyType = "ENSEMBL",
                            ont = "MF",
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 1,
                            qvalueCutoff = 1,
                            readable = T)
gene.GO_MF

gene.GO_CC <- enrichGO(gene = data$ID,
                            OrgDb = org.Dr.eg.db,
                            keyType = "ENSEMBL",
                            ont = "CC",
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 1,
                            qvalueCutoff = 1,
                            readable = T)
gene.GO_CC
####结果输出
gene.GO_CC_data <- data.frame(gene.GO_CC)
gene.GO_MF_data <- data.frame(gene.GO_MF)
gene.GO_BP_data <- data.frame(gene.GO_BP)
gene.GO_BP_data$type <- "BP"
gene.GO_MF_data$type <- "MF"
gene.GO_CC_data$type <- "CC"
gene.GO_ALL<-rbind(gene.GO_BP_data, gene.GO_MF_data)  ##### 列名一样直接合并用
gene.GO_ALL<-rbind(gene.GO_ALL, gene.GO_CC_data)#####  一次只能合并两个data.frame
write.table(gene.GO_ALL,"gene.GO_ALL_result",quote=F,sep = "\t")
###### kegg
search_kegg_organism("zebrafish", by="common_name")

gene.KEGG_low <- enrichKEGG(gene = data_low$V2,
                            organism = "dre",
                            keyType = "ncbi-geneid",
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
gene.KEGG_low

gene.KEGG <- enrichKEGG(gene = data$V2,
                            organism = "dre",
                            keyType = "ncbi-geneid",
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
gene.KEGG
#####作图
#go
##用lowlow的自带程序画???
high_BP<-dotplot(gene.GO_high_BP)
high_BP
high_MF <-dotplot(gene.GO_high_MF, showCategory = 10)
high_MF
high_CC <-dotplot(gene.GO_high_CC, showCategory = 10)
high_CC
plot_grid(high_BP,high_MF,high_CC,nrow = 3)
#######用ggplot画图
###将数据转化为data.frame,并按FDR排序，取前几???
#???
gene.GO_high_BP_head<-head(arrange(data.frame(gene.GO_high_BP),p.adjust),n=10)  
gene.GO_high_MF_head<-head(arrange(data.frame(gene.GO_high_MF),p.adjust),n=10) 
gene.GO_high_CC_head<-head(arrange(data.frame(gene.GO_high_CC),p.adjust),n=10) 
#???
gene.GO_BP_head<-head(arrange(data.frame(gene.GO_BP),p.adjust),n=10)  
gene.GO_MF_head<-head(arrange(data.frame(gene.GO_MF),p.adjust),n=10) 
gene.GO_CC_head<-head(arrange(data.frame(gene.GO_CC),p.adjust),n=10) 
 ####分数转小???
gene.GO_high_BP_head$GeneRatio <- parse_ratio(gene.GO_high_BP_head$GeneRatio)
gene.GO_high_MF_head$GeneRatio <- parse_ratio(gene.GO_high_MF_head$GeneRatio)
gene.GO_high_CC_head$GeneRatio <- parse_ratio(gene.GO_high_CC_head$GeneRatio)
#all
gene.GO_BP_head$GeneRatio <- parse_ratio(gene.GO_BP_head$GeneRatio)
gene.GO_MF_head$GeneRatio <- parse_ratio(gene.GO_MF_head$GeneRatio)
gene.GO_CC_head$GeneRatio <- parse_ratio(gene.GO_CC_head$GeneRatio)
#
gene.GO_high_BP_head$type <- "BP"
gene.GO_high_MF_head$type <- "MF"
gene.GO_high_CC_head$type <- "CC"
#all
gene.GO_BP_head$type <- "BP"
gene.GO_MF_head$type <- "MF"
gene.GO_CC_head$type <- "CC"
#
GO_high<-rbind(gene.GO_high_BP_head, gene.GO_high_MF_head)  ##### 列名一样直接合并用
GO_high<-rbind(GO_high, gene.GO_high_CC_head)#####  一次只能合并两个data.frame
GO_high$type=factor(GO_high$type, levels=c("CC","MF","BP"))
GO_term_order=factor(as.integer(rownames(GO_high)),labels=GO_high$Description)
#all
GO_all<-rbind(gene.GO_BP_head, gene.GO_MF_head)  ##### 列名一样直接合并用
GO_all<-rbind(GO_all, gene.GO_CC_head)#####  一次只能合并两个data.frame
GO_all$type=factor(GO_all$type, levels=c("CC","MF","BP"))
GO_term_order_all=factor(as.integer(rownames(GO_all)),labels=GO_all$Description)
#plot
COLS<-c("#66C3A5", "#8DA1CB", "#FD8D62")
ggplot(data=GO_high, aes(x=Count,y=GO_term_order, fill=type)) +
  geom_bar(stat="identity", width=0.8)+
  xlab("Num of Genes") + 
  ylab("GO term") + 
  scale_fill_manual(values=COLS)+
  theme_bw()

ggplot(data=GO_high,aes(x=100*(GeneRatio),y=GO_term_order))+ 
  geom_point(aes(size=Count))+
  geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
  scale_color_gradient(low="blue",high = "red")+
  labs(x="GeneRatio %",y="GO term",
       size="Num of Genes",color="-log10(FDR)")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())
plot_grid(p3,p4,ncol = 2)  
#all
p1<-ggplot(data=GO_all, aes(x=Count,y=GO_term_order_all, fill=type)) +
  geom_bar(stat="identity", width=0.8)+
  xlab("Num of Genes") + 
  ylab("GO term") + 
  scale_fill_manual(values=COLS)+
  theme_bw()


p2<-ggplot(data=GO_all,aes(-1*log10(GeneRatio),y=GO_term_order_all))+ 
  geom_point(aes(size=Count))+
  geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
  scale_color_gradient(low="blue",high = "red")+
  labs(x="-log10 GeneRatio",y="GO term",
  size="Num of Genes",color="-log10(FDR)")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())
plot_grid(p1,p2,ncol = 2)  


cnetplot(gene.GO_high, showCategory = 5)
emapplot(gene.GO_high, showCategory = 20)



heatplot(gene.GO_high, foldChange =data_high$V5)
###kegg
barplot(gene.KEGG_low,title="Tnf-ab low experss KEGG Enrichment")


ggplot(data=gene.KEGG_low, aes(x=Count,y=Description, fill=p.adjust)) +
  geom_bar(stat="identity")+
  xlab("Num of Genes") + 
  ylab("KEGG term")


########## GSEA  enrichment
#预处??? 排序
setwd("C:\\Users\\lcc\\Desktop\\ear")
geneList <- select(data2,ID,SYMBOL,logFC,FDR)
geneList.sort <- arrange(geneList,desc(logFC),FDR); head(geneList.sort)

geneList.sort$SYMBOL<-toupper(geneList.sort$SYMBOL)##基因名换为大???
gene <- geneList.sort$SYMBOL
####读取broad GSEA提供的gene sets
c5 <- read.gmt("c5.all.v7.1.symbols.gmt");head(c5)
h <-read.gmt("h.all.v7.1.symbols.gmt");head(h)
c2 <-read.gmt("c2.all.v7.1.symbols.gmt");head(c2)
enrich_c5 <- enricher(gene, TERM2GENE=c5); head(enrich_c5)
enrich_h <- enricher(gene, TERM2GENE=h); head(enrich_h)
enrich_c2 <- enricher(gene, TERM2GENE=c2); head(enrich_c2)

glist <- geneList[,3];head(glist)## feature 1: numeric vector
names(glist) <- as.character(geneList[,2]);head(glist)## feature 2: named vector
glist <- sort(glist,decreasing = T); head(glist)## feature 3: decreasing order

gsea <- GSEA(glist, TERM2GENE=c5,pvalueCutoff =1); head(gsea)

gsea.go <- gseGO(glist,OrgDb = org.Dr.eg.db, pvalueCutoff =1); head(gsea.go)
gseaplot2(gsea, 1)
####2

list<-c("GO:0007250","GO:0007253","GO:0032088","GO:0033256","GO:0033257","GO:0035525","GO:0038061","GO:0038167","GO:0038168","GO:0039644","GO:0039652","GO:0043122","GO:0043123","GO:0043124","GO:0051059","GO:0051092","GO:0061765","GO:0071159","GO:0085032","GO:0085033","GO:0085034","GO:1901222","GO:1901223","GO:1901224")


geneList2 <- select(data2,ID,logFC)
geneList2.sort <- arrange(geneList2,desc(logFC)); head(geneList.sort)

write.table(GO_DATA$PATHID2NAME,"GO_DATA_PATHID2NAME",quote=F,sep = "\t")




GO_0007249<-unlist(GO_DATA$PATHID2EXTID["GO:0007249"])
GO_0007249<-cbind(rep("NF-kappaB",length(GO_0007249)),as.data.frame(GO_0007249))
colnames(GO_0007249)<-c("ont","gene")

GO_0007250<-unlist(GO_DATA$PATHID2EXTID["GO:0007250"])
GO_0007250<-cbind(rep("NF-kappaB",length(GO_0007250)),as.data.frame(GO_0007250))
colnames(GO_0007250)<-c("ont","gene")

NFKB<- rbind(GO_0007249,GO_0007250)



for (i in 1:length(list)){
  go <-list[i]
  #name<-past("go",i,sep="")
  NFKB_i<-unlist(GO_DATA$PATHID2EXTID["GO:0007250"])
  NFKB_i<-cbind(rep("NF-kappaB",length(NFKB_i)),as.data.frame(NFKB_i))
  colnames(NFKB_i)<-c("ont","gene")
  NFKB<- rbind(NFKB,NFKB_i)

}


head(NFKB)

gsea <- GSEA(glist, TERM2GENE=GO_DATA, verbose=FALSE,pvalueCutoff =1); head(gsea)
gseaplot(gsea, 1)
gseaplot2(gsea, 1)
