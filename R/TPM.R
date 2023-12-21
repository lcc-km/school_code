library(ggplot2)
library(grid)

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

####################阳性位点的基因##############################################
list_goal<-c('ENSDARG00000024433','ENSDARG00000099517')
setwd("D:\\backup\\our_data\\RNA\\TPM")

########################包含又raw count数据的RNA-seq结果################
list.files<-c('RMA-1_FRRB190272334')   #,'RMA-2_FRRB190272335','RMA-3_FRRB190272336','RMA-4_FRRB190272337',
#'RMA-5_FRRB190272338','RMA-6_FRRB190272339','RMB-1_FRRB190272340','RMB-2_FRRB190272341',
#'RMB-3_FRRB190272342','RMB-4_FRRB190272343','RMB-5_FRRB190272344','RMB-6_FRRB190272345',
#'RMC-1_FRRB190272346','RMC-2_FRRB190272347','RMC-3_FRRB190272348','RMC-4_FRRB190272349',
#'RMC-5_FRRB190272350','RMC-6_FRRB190272351','RMD-1_FRRB190272352','RMD-2_FRRB190272353',
#'RMD-3_FRRB190272354','RMD-4_FRRB190272355','RMD-5_FRRB190272356','RMD-6_FRRB190272357',
#'RME-1-2_FRRB190276052','RME-2-2_FRRB190276053','RME-3_FRRB190272360','RME-4_FRRB190272361',
#'RME-5_FRRB190272362','RME-6_FRRB190272363','RMF-1_FRRB190272364','RMF-2_FRRB190272365',
#'RMF-3_FRRB190272366','RMF-4_FRRB190272367','RMF-5_FRRB190272368','RMF-6_FRRB190272369'
for (i in 1:length(list.files)) {
read_file_name  <-paste(list.files[[i]],".txt",sep="")
countdata<-read.table(read_file_name ,sep="\t",header = FALSE)
countDf <- data.frame(id=countdata$V1,count=countdata$V7,length=countdata$V6)  #######选出需要的�?

countDf$effLength <- countDf$length - 150 +1    ##########不懂值怎么来的

countDf <- countDf[countDf$effLength>1,]

####### 为刚刚的公式赋�?
countDf$tpm <- with(countDf, countToTpm(count, effLength))
countDf$fpkm <- with(countDf, countToFpkm(count, effLength))
with(countDf, all.equal(tpm, fpkmToTpm(fpkm)))
countDf$effCounts <- with(countDf, countToEffCounts(count, length, effLength))

#########选取阳性基�?
count_gola <- countDf[countDf$id%in%list_goal,]
count_gola_work <- data.frame(count_gola$id,count_gola$tpm)
######### 计算TPM平均�?
mean<-mean(countDf$tpm)
##########  绘图
ggplot(count_gola_work,aes(x=count_gola.id,y=count_gola.tpm))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = mean)

p_split1 <- p+ coord_cartesian(ylim = c(0, 100))+theme(axis.title.y = element_blank())

p_split2 <- p+ coord_cartesian(ylim = c(150, 500))+theme(axis.text.x = element_blank(),axis.title.y = element_blank(),
                                                         axis.ticks.x = element_blank(),axis.title.x=element_blank())
p_split3 <- p+ coord_cartesian(ylim = c(10000, 30000))+theme(axis.text.x = element_blank(),axis.title.y = element_blank(),
                                                             axis.ticks.x = element_blank(),axis.title.x=element_blank())


p_name <- paste(list.files[[i]],".pdf",sep="")


#pdf(p_name)
#grid.newpage()
plot_site1 <- viewport(x = 0.02, y = 0, width = 1, height = 0.5, just = c('left', 'bottom'))
plot_site2 <- viewport(x = 0.02, y = 0.5, width = 1, height = 0.2, just = c('left', 'bottom'))
plot_site3 <- viewport(x = 0, y = 0.7, width = 1.02, height = 0.2, just = c('left', 'bottom'))
print(p_split1, vp = plot_site1)
#print(p_split2, vp = plot_site2)
#print(p_split3, vp = plot_site3)
#dev.off()
}

