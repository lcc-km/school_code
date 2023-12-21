## This pipline is used to compare gene expression among multiple samples,
## including: 1)one interal control gene  2) one reference sample
require("ggplot2")
require("ggsignif")
require("ggpubr")
library(agricolae)
#################################################################################
########## !! Define input informations below: ##################################

## prepare a ""TAB seprated" file has at least 3 colums: "gene", "sample" ,"ct"
## change the following file directory & file names accordingly
setwd("D:\\backup\\enhancer paper\\qPCR\\wj_text")  # file directory
d <- read.csv("qPCR.txt",sep="\t",header = T) # file name

# check name of columns, which should include : "gene", "sample" ,"ct". 
# If NOT, change the input file.
names(d)  


## setup names of genes and samples:
levels(factor(d$gene)) 
levels(factor(d$sample))
ref_gene <- "elf"     # give name of refrence gene (internal control gene)
ref_sample <- "36hpf"  # give name of CONTROL sample
# fill <- c("gray60","gray15") #set color

# set order of genes:
gene_order <- c("tgfa","tgfb","tgfbr2b",
                "tp53","cdkn1a","cdkn1bb","rb1","il6","il1b")

# set order of samples:
sample_order <- c("36hpf", "3dpf" , "5dpf" )

# check SD, repeat experiment if sd >0.5
d$group <- paste(d$gene,d$sample,sep ="_" )
sd <- tapply(d$ct,INDEX = d$group,FUN = sd)
sort(sd,decreasing = T)   ##  check standard diviation




###############################################################################
## !! DO NOT CHANGE from here except for advanced users:
##############################################################################

# 1. data preparation
  ## internal mean: 
  dref <- data.frame(ref_ct=tapply(d[d$gene==ref_gene,]$ct,
                                     INDEX = d[d$gene==ref_gene,]$group,FUN = mean))
  dref$sample <- sapply(strsplit(as.character(row.names(dref)),'_'), "[", 2) ## seprated to get group name
  d1 <- merge(d,dref,by="sample")

  ## (minus) delta ct:  reference-target ######
  d1$delta_ct <- d1$ref_ct-d1$ct
  d1 <- d1[d1$gene != ref_gene,]

# 2. statistic of one-way ANOVA and multiple comparison: TukeyHSD within gene: 
  a <- list()
  HSD_w_r <- data.frame(delta_ct=NA,groups=NA,sample=NA,gene=NA)
  for(g in gene_order){
    aov.d <- d1[d1$gene==g,]
    fit <- aov(delta_ct~sample,aov.d)
    HSD<-HSD.test(fit,"sample")
    HSD_w<-HSD$groups
    HSD_w$sample<-row.names(HSD_w)
    HSD_w$gene <- g
    HSD_w_r <- rbind(HSD_w_r,HSD_w)
    hsd <- TukeyHSD(fit)
    names(hsd) <- g
    a <- list(title="Results of TukeyHSD multiple comparison within gene",
              result=c(a$result,hsd))
  }
  a # show results of TukeyHSD multiple comparison
  HSD_w_r
  HSD_w_r <- na.omit(HSD_w_r)
  HSD_w_r$work <-paste(HSD_w_r$sample,HSD_w_r$gene,sep = "_")
  ## average delta CT & data for plot:
  mean <- tapply(d1$delta_ct,INDEX = d1$group,FUN = mean)
  sem <- function(x) sd(x)/sqrt(length(x)) # define function for standard error of mean
  se <- tapply(d1$delta_ct,INDEX = d1$group,FUN = sem)
  
  dp <- data.frame(mean,se)
  dp$group <- row.names(dp)
  dp$gene <- sapply(strsplit(as.character(row.names(dp)),'_'), "[", 1)
  dp$sample <- sapply(strsplit(as.character(row.names(dp)),'_'), "[", 2)
  
  
# 3.delta delta Ct plot normlized to ref sample  (No Significance marker)
dp1 <- dp[!(dp$gene==ref_gene),]
dp1$gene <-factor(dp1$gene,levels = gene_order)
dp1$sample <- factor(dp1$sample,levels = sample_order)

##### delta delta ct:  delta ct - ref sample ######
dp1$delta_delta_ct <- 0
for (i in levels(dp1$gene)) {
  ctrl <- mean(as.numeric(dp1[dp1$sample==ref_sample & dp1$gene==i,"mean"]))
  dp1[dp1$gene==i,]$delta_delta_ct <- dp1[dp1$gene==i,]$mean-ctrl
}

dp1$norm_delta <- 2^dp1$delta_delta_ct
dp1$work <- paste(dp1$sample,dp1$gene,sep = "_")

dp2 <- merge(dp1,HSD_w_r,by="work")
  
ggplot(dp2,aes(x=sample.x,fill=sample.x))+
  geom_bar(aes(y=norm_delta),stat = "identity",width=0.8,color="black")+
  geom_errorbar(aes(ymax=norm_delta+2^se-1,ymin=norm_delta),width=0.6)+
  ylab(label = "2^ delta delta Ct to ref sample")+
  scale_fill_grey(start=0.0,end=0.99)+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.text.x=element_text(angle=90))+  # change x axis angle
  facet_wrap(~gene.x,nrow = 1)+
  geom_text(label=dp2$groups,y=((dp2$norm_delta+0.1)*2^dp1$se))+
  xlab("sample")+labs(fill="sample") 
  
  


