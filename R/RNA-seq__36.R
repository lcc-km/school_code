library(Vennerable)
setwd("C:\\Users\\lcc\\Desktop\\our_data\\RNA")
BA <- read.csv("zebBiA_exact_FC.txt",header=TRUE,sep = "\t")
BA_elect <- subset(BA,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

CA <- read.csv("zebCiA_exact_FC.txt",header=TRUE,sep = "\t")
CA_elect <- subset(CA,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

CB <- read.csv("zebCiB_exact_FC.txt",header=TRUE,sep = "\t")
CB_elect<- subset(CB,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

DA <- read.csv("zebDiA_exact_FC.txt",header=TRUE,sep = "\t")
DA_elect<- subset(DA,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

DB <- read.csv("zebDiB_exact_FC.txt",header=TRUE,sep = "\t")
DB_elect<- subset(DB,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

DC <- read.csv("zebDiC_exact_FC.txt",header=TRUE,sep = "\t")
DC_elect<- subset(DC,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

EA <- read.csv("zebEiA_exact_FC.txt",header=TRUE,sep = "\t")
EA_elect<- subset(EA,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

EB <- read.csv("zebEiB_exact_FC.txt",header=TRUE,sep = "\t")
EB_elect<- subset(EB,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

EC <- read.csv("zebEiC_exact_FC.txt",header=TRUE,sep = "\t")
EC_elect<- subset(EC,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

ED <- read.csv("zebEiD_exact_FC.txt",header=TRUE,sep = "\t")
ED_elect<- subset(ED ,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

FA <- read.csv("zebFiA_exact_FC.txt",header=TRUE,sep = "\t")
FA_elect<- subset(FA,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

FB <- read.csv("zebFiB_exact_FC.txt",header=TRUE,sep = "\t")
FB_elect<- subset(FB,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

FC <- read.csv("zebFiC_exact_FC.txt",header=TRUE,sep = "\t")
FC_elect<- subset(FC,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

FD <- read.csv("zebFiD_exact_FC.txt",header=TRUE,sep = "\t")
FD_elect<- subset(FD,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

FE <- read.csv("zebFiE_exact_FC.txt",header=TRUE,sep = "\t")
FE_elect<- subset(FE,abs(logFC)>=1,select = c(ID,logFC,logCPM,PValue,FDR))

all <- Reduce(intersect,list(v1=BA_elect$ID,v2=CA_elect$ID,v3=CB_elect$ID,v4=DA_elect$ID,
                      v5=DB_elect$ID,v5=DC_elect$ID,v7=EA_elect$ID,v8=EB_elect$ID,
                      v9=EC_elect$ID,v10=ED_elect$ID,v11=FA_elect$ID,v12=FB_elect$ID,
                      v13=FC_elect$ID,v14=FD_elect$ID,v15=FE_elect$ID))
all1 <- list(BA_elect$ID,CA_elect$ID,CB_elect$ID,DA_elect$ID,DB_elect$ID,
             DC_elect$ID,EA_elect$ID,EB_elect$ID,EC_elect$ID,ED_elect$ID,
             FA_elect$ID,FB_elect$ID,FC_elect$ID,FD_elect$ID,FE_elect$ID)

BA_ID <- data.frame(BA_elect$ID)
CA_ID <- data.frame(CA_elect$ID)
CB_ID <- data.frame(CB_elect$ID)
DA_ID <- data.frame()


