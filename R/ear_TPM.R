##############复制来的算法 #####################
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}
#########读取文件
setwd("C:\\Users\\lcc\\Desktop\\ear")
list.files<-c("AB-1-AB-2_FRRB20H100363-1a.txt_TPM",
              "AB-1_FRRB202333500-1a.txt_TPM",
              "AB-2_FRRB202333501-1a.txt_TPM",
              "Tnf-1_FRRB202333502-1a.txt_TPM",
              "Tnf-1-Tnf-2_FRRB20H100364-1a.txt_TPM",
              "Tnf-2_FRRB202333503-1a.txt_TPM")
for (i in 1:length(list.files)){
read_file_name  <-paste(list.files[[i]],sep="")
countdata<-read.table(read_file_name ,sep="\t",header =TRUE)
countDf <- data.frame(id=countdata$Geneid,count=countdata$Count,length=countdata$Length)  #######选出需要的行

countDf$effLength <- countDf$length - 150 +1    ##########不懂值怎么来的

countDf <- countDf[countDf$effLength>1,]

####### 为刚刚的公式赋值
countDf$tpm <- with(countDf, countToTpm(count, effLength))
countDf$fpkm <- with(countDf, countToFpkm(count, effLength))
with(countDf, all.equal(tpm, fpkmToTpm(fpkm)))
#countDf$effCounts <- with(countDf, countToEffCounts(count, length, effLength))
############结果输出
file_name <- paste(list.files[[i]],"TPM.result.txt",sep="")
write.table(countDf,file_name,quote=F,sep = "\t")
}
