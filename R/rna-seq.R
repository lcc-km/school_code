setwd("C:\\Users\\49248\\Desktop\\parallel comparisonRNA-seq\\28-zeb\\muscle")
a <- read.csv("zebmuscle-brain_exact_FC1_FDR0.05.txt",header=TRUE,sep = "\t")
a_elect <- subset(a,logFC>=1,select = c(ID,logFC,logCPM,PValue,FDR))
b <- read.csv("zebmuscle-gill_exact_FC1_FDR0.05.txt",header=TRUE,sep = "\t")
b_elect <- subset(b,logFC>=1,select = c(ID,logFC,logCPM,PValue,FDR))
c <- read.csv("zebmuscle-heart_exact_FC1_FDR0.05.txt",header=TRUE,sep = "\t")
c_elect<- subset(c,logFC>=1,select = c(ID,logFC,logCPM,PValue,FDR))
d <- read.csv("zebmuscle-intestine_exact_FC1_FDR0.05.txt",header=TRUE,sep = "\t")
d_elect<- subset(d,logFC>=1,select = c(ID,logFC,logCPM,PValue,FDR))
e <- read.csv("zebmuscle-kidney_exact_FC1_FDR0.05.txt",header=TRUE,sep = "\t")
e_elect<- subset(e,logFC>=1,select = c(ID,logFC,logCPM,PValue,FDR))
f <- read.csv("zebmuscle-liver_exact_FC1_FDR0.05.txt",header=TRUE,sep = "\t")
f_elect<- subset(f,logFC>=1,select = c(ID,logFC,logCPM,PValue,FDR))

x <- Reduce(intersect,list(v1=a_elect$ID,v2=b_elect$ID,v4=d_elect$ID,v6=f_elect$ID))
y <- Reduce(intersect,list(v1=a_elect$ID,v2=b_elect$ID,v3=c_elect$ID,v4=d_elect$ID,v5=e_elect$ID,v6=f_elect$ID))
z <- Reduce(intersect,list(v1=a_elect$ID,v2=b_elect$ID,v3=c_elect$ID,v4=d_elect$ID,v6=f_elect$ID))
filename1 <- paste(x)
filename2 <- paste(y)
filename3 <- paste(z)
write.csv(filename1,"C:\\Users\\49248\\Desktop\\parallel comparisonRNA-seq\\28-zeb\\muscle\\result_normal_no_brain_Kidney.csv",quote = F)
write.csv(filename2,"C:\\Users\\49248\\Desktop\\parallel comparisonRNA-seq\\28-zeb\\muscle\\result_normal.csv",quote = F)
write.csv(filename3,"C:\\Users\\49248\\Desktop\\parallel comparisonRNA-seq\\28-zeb\\muscle\\result_normal_no_brain.csv",quote = F)
