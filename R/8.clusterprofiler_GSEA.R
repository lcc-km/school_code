library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(KEGG.db)
library(enrichplot)
library(dplyr)


#####################################################################################################
BiocManager::install("org.Mm.eg.db")

####   DO NOT CHANGE FROM HERE??
setwd("C:\\Users\\admin\\OneDrive - st.shou.edu.cn\\桌面")
dd <- read.table("cochlearE14-P33_exact_FC1_FDR0.05.txt",header=T,sep = "\t")
dd$gene<-rownames(dd)



keytypes(org.Mm.eg.db)
ENTREZID <- bitr(dd$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
ENTREZID$ENTREZID
colnames(ENTREZID)<-c("gene","ENTREZID")
#colnames(all_id)=c("gene","NCBIID","emD")

d2 <- merge(dd,ENTREZID,by.x="gene")
d2<- arrange(d2,desc(logFC))

## 1. GSEA: KEGG  #######################################################
  # rank genelist
  d2$fcsign <- sign(d2$logFC)
  d2$logP <- -log10(d2$PValue)
  d2$metric <-  d2$logP/d2$fcsign
  GSEA_input<-d2$metric
  names(GSEA_input) = as.character(d2$ENTREZID)
  GSEA_input = sort(GSEA_input, decreasing = T)

gseKEGG.res <- gseKEGG(GSEA_input, organism = "mmu",keyType = "ncbi-geneid",nPerm = 1000, minGSSize = 5, maxGSSize = 1000, pvalueCutoff=1)
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2\\sc_RNA_edgeR\\gsea")
write.table(as.data.frame(gseKEGG.res),file=paste(prefix,"GSEA.kegg.enrich.ncbi.txt",sep = "_"),quote = F,sep="\t",row.names = F)


inpath1 <-"mmu04310"	#Wnt signaling pathway

inpath2<-"mmu04064" #	NF-kappa B signaling pathway   
inpath3<-"mmu04668" #	TNF signaling pathway
inpath4<-"mmu04115" # 	p53 signaling pathway
inpath5<-"mmu04010" #	MAPK signaling pathway
inpath6<-"mmu04330" #	Notch signaling pathway
inpath7<-"mmu04150" #mTOR signaling pathway
inpath8<-"mmu04151" #	PI3K-Akt signaling pathway


inpath
p1 <- gseaplot2(gseKEGG.res,geneSetID = inpath1,pvalue_table = T)
p1
ggsave(paste(prefix,inpath,"GSEA.pdf",sep="_"),p1,width = 6,height = 5)  

core <- strsplit(gseKEGG.res@result[gseKEGG.res@result$ID==inpath,]$core_enrichment,split = "/")[[1]]
bitr(geneID = core,fromType = "ENTREZID",toType ="SYMBOL" ,OrgDb = org.Dr.eg.db)

##### GSEA: GO

gseGO.res <- gseGO(GSEA_input,ont = "BP",OrgDb = org.Dr.eg.db,keyType = "ENTREZID",pvalueCutoff = 1 )
write.table(as.data.frame(gseGO.res),file=paste(prefix,"GSEA.GOBP.enrich.ncbi.txt",sep = "_"),quote = F,sep="\t",row.names = F)

inpath <-"GO:0090596" #	sensory organ morphogenesis
inpath <-"GO:0050896" #	response to stimulus
inpath <-"GO:0023052"	#signaling
inpath <-"GO:0030154"	#cell differentiation
inpath <-"GO:0072358" #	cardiovascular system development
inpath <-"GO:0060326" #	cell chemotaxis
inpath <-"GO:0009605" #	response to external stimulus
p2 <- gseaplot2(gseGO.res,geneSetID = inpath,pvalue_table = T)
p2
ggsave(paste(prefix,strsplit(inpath,":")[[1]][2],"GSEA.GOBP.pdf",sep="_"),p2,width = 6,height = 5)

core <- strsplit(gseGO.res@result[gseGO.res@result$ID==inpath,]$core_enrichment,split = "/")[[1]]
bitr(geneID = core,fromType = "ENTREZID",toType ="SYMBOL" ,OrgDb = org.Dr.eg.db)


####  GSEA ??? GO CC
gseGOCC.res <- gseGO(GSEA_input,ont = "CC",OrgDb = org.Dr.eg.db,keyType = "ENTREZID",pvalueCutoff = 1 )
write.table(as.data.frame(gseGOCC.res),file=paste(prefix,"GSEA.GOCC.enrich.ncbi.txt",sep = "_"),quote = F,sep="\t",row.names = F)

