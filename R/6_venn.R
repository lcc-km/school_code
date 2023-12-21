library(venn)
setwd("D:\\backup\\our_data\\RNA\\raw_clean\\vcf\\WORK\\gene.locate\\F")
F <- read.csv("F.result",header=FALSE,sep = "\t")
F1 <- read.csv("RMF-1_FRRB190272364.raw.vcf_result.txt",header=FALSE,sep = "\t")
F2 <- read.csv("RMF-2_FRRB190272365.raw.vcf_result.txt",header=FALSE,sep = "\t")
F3 <- read.csv("RMF-3_FRRB190272366.raw.vcf_result.txt",header=FALSE,sep = "\t")
F4 <- read.csv("RMF-4_FRRB190272367.raw.vcf_result.txt",header=FALSE,sep = "\t")
F5 <- read.csv("RMF-5_FRRB190272368.raw.vcf_result.txt",header=FALSE,sep = "\t")
F6 <- read.csv("RMF-6_FRRB190272369.raw.vcf_result.txt",header=FALSE,sep = "\t")

FF <- paste(F$V3,F$V4,F$V5,F$V6,sep="_")
F1_locate <- paste(F1$V3,F1$V4,F1$V5,F1$V6,sep="_")
F2_locate <- paste(F2$V3,F2$V4,F2$V5,F2$V6,sep="_")
F3_locate <- paste(F3$V3,F3$V4,F3$V5,F3$V6,sep="_")
F4_locate <- paste(F4$V3,F4$V4,F4$V5,F4$V6,sep="_")
F5_locate <- paste(F5$V3,F5$V4,F5$V5,F5$V6,sep="_")
F6_locate <- paste(F6$V3,F6$V4,F6$V5,F6$V6,sep="_")
Funique <- unique(FF)
write.table(Funique,"Funique",sep = "\t")
data <- list(F1=F1_locate,F2=F2_locate,F3=F3_locate,F4=F4_locate,F5=F5_locate,F6=F6_locate)

venn(data,zcolor='style')


