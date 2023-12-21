library(ggplot2)
library(qvalue)
setwd("D:\\backup\\our_data\\RNA\\VQSR\\10\\gene_area\\expression_quantity\\integrate_locate\\all\\logMm\\ID\\wj")


file.name <- 'ENSDARG00000017441'
limit <- 10000000  # 文件中存在非预期范围内的位点，剔除该位点使用的标准
p.thrd <- 0.05  # 统计检验的显著标准

## 读入 fish的genotype 预测结果
d <- read.csv(file = file.name,sep = "\t",header = F)
site <- unique(d$V16)

## 获取 RNA的reads 并计算log AEI
a <- unique(d[,c(2,6,7,8,9)])
colnames(a) <- c("id","a1","a2","an1","an2")
a$log_AEI <- abs(log2(a$an1/a$an2))
a$fish_id <- substr(a$id,1,2)

a1 <- a[,c("fish_id","log_AEI")]

## 选择genome上一个SNP位点 进行测试
s1 <- d[d$V16==site[110],c(2,4,16,19,17,22,20,24)]  ## 选择其中一个SNP位点
s1
colnames(s1) <- c("id","chr","pos","p1","allele1","p2","allele2","genotype")
s1$geno2 <- ifelse(s1$allele1==s1$allele2,0,1)  ## 再次计算genotype
s1

s1$fish_id <- substr(s1$id,1,2)
d.plot <- merge(a1,s1,by = "fish_id") # 绘图用数据

x <- d.plot$genotype
y <- d.plot$log_AEI

test <- wilcox.test(y~x,alternative="less") # 使用单尾检验，Ha： group 0 < group 1
test
test$p.value

p1<-ggplot(d.plot,aes(x=as.factor(genotype),y=log_AEI))+
  geom_boxplot()+
  geom_jitter(width = 0.2)
name1_1<- paste(file.name,"_p1",".png",sep="")
ggsave(file=name1_1,plot=p1)


#####   单个位点 end #########


## 生成所有位点的phenotype 及 genotype 数据：SS为最终的理想数据格式,map为snp位置数据
ss <- a1
map <- data.frame(id=NA,chr=NA,pos=NA)
map<- map[-1,]

i=1
for (i in 1:length(site)) {
  s1 <- d[d$V16==site[i],c(2,4,16,19,17,22,20,24)]  ## 选择其中一个SNP位点
  s1
  colnames(s1) <- c("id","chr","pos","p1","allele1","p2","allele2","genotype")
  s1$geno2 <- ifelse(s1$allele1==s1$allele2,0,1)  ## 再次计算genotype
  s1

  df <- data.frame(id=paste("s",i,sep=""),chr=unique(s1$chr),pos=unique(s1$pos))
  map <- rbind(map,df)
}
map # SNP 位置信息文件，类map格式文件
write.table(map,file = paste(file.name,".map","_2",sep=""),quote = F,sep="\t",row.names = F)
gc()

print(length(site))
i=1
for (i in 1:length(site)) {
  s1 <- d[d$V16==site[i],c(2,4,16,19,17,22,20,24)]  ## 选择其中一个SNP位点
  s1
  colnames(s1) <- c("id","chr","pos","p1","allele1","p2","allele2","genotype")
  s1$geno2 <- ifelse(s1$allele1==s1$allele2,0,1)
  #name1 <- paste("geno",unique(s1$chr),unique(s1$pos),sep = "_")
  name1 <- paste("s",i,sep="")
  colnames(s1)[length(colnames(s1))] <-name1

  s1$fish_id <- substr(s1$id,1,2)
  ns1 <- ncol(s1)

  ss <- merge(ss,s1[,c(ns1,ns1-1)],by="fish_id",all = T)
  ss=ss[!duplicated(ss),]
  print(i)


}
ss  # SNP 的genotype 和 phenotype文件，类似ped格式文件

write.table(ss,file = paste(file.name,".ped.het1","_2",sep=""),quote = F,sep="\t",row.names = F)

# 统计检验，获得所有位点的检验p值，并绘图（此例使用Wilcoxon rank sum test，对满足正态的数据，可是使用t test）
## 生成一个空dataframe
t.d <- data.frame(chr=NA,pos=NA,t.pval=NA,w.pval=NA)
t.d<- t.d[-1,]

for (i in 1:(ncol(ss)-2)) {
  d.p <- ss[,c(1,2,i+2)]
  geo.name <-colnames(d.p)[3]
  # chr <- strsplit(geo.name,split = "_")[[1]][2] ## 获取染色体号
  # pos <- strsplit(geo.name,split = "_")[[1]][3] ## 获取SNP位置号
  chr <- map[map$id==geo.name,"chr"] ## 获取染色体号
  pos <- map[map$id==geo.name,"pos"] ## 获取SNP位置号

  colnames(d.p)[3] <- "genotype"
  x <- d.p$genotype
  y <- d.p$log_AEI

  level <- length(unique(na.omit(x)))

  if (level<2) {  ## 对于一些 去除NA后只剩一种genotype的位点，需要去除
    tp <- NA
    wp <- NA

  }else {
    tt <- t.test(y~x,alternative="less") ## 使用单尾检验，Ha： group 0 < group 1
    tp <- tt$p.value

    wt <- wilcox.test(y~x,alternative="less") # 使用单尾检验，Ha： group 0 < group 1
    wp <- wt$p.value
  }


  #
  df <- data.frame(chr=chr,pos=pos,t.pval=tp,w.pval=wp)
  t.d <- rbind(t.d,df)

  d.p1 <- na.omit(d.p) ## 去除带有na的行
  # ggplot(d.p1,aes(x=as.factor(genotype),y=log_AEI))+
  #   geom_boxplot()+
  #   geom_jitter(width = 0.2)

}
t.d
t.d <- na.omit(t.d)
t.d <- t.d[t.d$pos>limit,]  ###  剔除位置不正确的位点

# t test 结果
t.d$t.fdr <-  p.adjust(t.d$t.pval,"BH") # global FDR矫正
t.d$t.qval <- qvalue(t.d$t.pval)$qvalue      # qvalue矫正

# wilcox rank sum test 结果：
t.d$w.fdr <-  p.adjust(t.d$w.pval,"BH") # global FDR矫正
t.d$w.qval <- qvalue(t.d$w.pval)$qvalue      # qvalue矫正

min(t.d$t.fdr,na.rm = T)
min(t.d$t.qval,na.rm = T)
##  t.d的是最主要的结果文件

#导出文件
write.table(t.d,file = paste(file.name,".p",sep=""),quote = F,sep="\t",row.names = F)


## 针对 log AEI绘图：

p2<-ggplot(ss[,1:2],aes(x=log_AEI))+
  geom_density(alpha=.2, fill="#FF6666")+
  labs(title=file.name)
name2<- paste(file.name,"_p2",".png",sep="")
ggsave(file=name2,plot=p2)




name3<- paste(file.name,"_p3",".png",sep="")
png(file=name3)
qqnorm(ss$log_AEI) # 绘qq图
dev.off()



#绘图使用wilcox fdr：(图可以使用ggsave 导出)
t.d$col <- "black"
#t.d[t.d$w.fdr<p.thrd,]$col <-"red"


mycolors<-c("black","red")
p4<-ggplot(t.d,aes(x=pos,y=-log10(t.d$w.fdr)))+
  geom_point(aes(color=as.factor(col)),shape=1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
  scale_color_manual(values = mycolors)+
  geom_hline(yintercept = -log10(p.thrd),colour="grey50")+
  ylab("-log10 FDR")+labs(title=file.name)
name4<- paste(file.name,"_p4",".png",sep="")
ggsave(file=name4,plot=p4)

print(file.name)



