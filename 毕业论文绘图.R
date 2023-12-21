library(ggplot2)
#BiocManager::install("rlang")
setwd("D:\\backup\\毕业论文\\")
# 含灰色的调色???
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#####RNA比对???
data_1 <-read.table("RNA_mapping.txt",header = T,sep = "\t")
data_1$Compare.rates_2 <- as.numeric(substr(data_1$Compare.rates,1,5))####提取数字并转化为数字类型
data_1$samples_name <-substr(data_1$samples.2,3,5)
data_1$class <-factor(substr(data_1$samples_name,1,1)) 
p<- ggplot(data_1,aes(samples_name,Compare.rates_2,fill=class))+
  geom_bar(stat = "identity",width=0.5,fill="#CC79A7")+   
  theme_classic()##去除背景
p+  ylim(0,100)  ####改变y轴范???
p+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  guides(fill=FALSE) +  ###旋转坐标轴label的方???,取消图例
  xlab("F1 individual")+
  ylab("Compare rates %")+
  theme(axis.text.x = element_text(size=20),axis.text.y=element_text(size=20),
        axis.title.x= element_text(size=15),axis.title.y= element_text(size=20))
   

########DNA比对???
data_2 <-read.table("DNA_mapping.txt",header = T,sep = "\t")
ggplot(data_2,aes(sample_2,rate,fill=class))+
  geom_bar(position = "dodge", stat = "identity",width=0.5)+  ###分组参数
  theme_classic()+ylim(0,100) +
  xlab("F0 individual")+ylab("Compare rates %")+
  theme(legend.position="bottom") +##改变图例的位???,左边left,右边 right, 底部bottom
  scale_fill_manual(values=c("#CC79A7","#56B4E9"))+  ###设置分列颜色
  theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=20),
        axis.title.x= element_text(size=20),axis.title.y= element_text(size=20))
  
####RNA-seq 具有相同基因相同位点的数???
#######基因组DP——qc
setwd("G:\\backup\\毕业论文")
data <- read.table("基因组_DP_QC.txt",header = T)
ggplot(data,aes(x=sample,y=number ,fill=QC))+
  geom_bar(position = "dodge", stat = "identity",width=0.5)+  ###分组参数
  theme_classic()+
  xlab("F0 individual")+ylab("Compare rates %")+
  theme(legend.position="bottom") +##改变图例的位???,左边left,右边 right, 底部bottom
  scale_fill_manual(values=c("#CC79A7","#56B4E9"))+
  geom_text(aes(label=number))
####DNA数据中的SNP与RNA-seq中一致的比例
library(grid)
setwd("G:\\backup\\our_data\\RNA\\VQSR\\10\\gene_area\\expression_quantity\\integrate_locate\\all\\logMm\\RNA_DNA_coherence")
data_3 <- read.csv("G:\\backup\\our_data\\RNA\\VQSR\\10\\gene_area\\expression_quantity\\integrate_locate\\all\\logMm\\RNA_DNA_coherence\\coherence_ratio.txt",sep = "\t",header = TRUE)
ggplot(data_3,aes(x=class,y=Genetic_ratio,fill=class))+
  geom_boxplot()+theme(legend.position="top")+ylim(80, 100)+theme_classic()+   ###去除背景网格
  ylab("Genetic ratio %")+guides(fill=FALSE)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7"))

####RNA-seq data count
data_4<- read.table("G:\\backup\\our_data\\RNA\\VQSR\\10\\gene_area\\expression_quantity\\RNA-seq_data_count_2.txt",header=TRUE,sep="\t")
data_4$type<-factor(data_4$type,levels=c("all_SNP","aftet_VQSR","heterozygosis","locate_in_gene"),ordered=TRUE)
ggplot(data_4,aes(x=class,y=number,width=1))+
  geom_boxplot(aes(fill=type))+theme(legend.position="top") +
  theme_classic()+
  scale_fill_grey()+
  theme(legend.position="bottom",axis.text.x = element_text(size=15),axis.text.y=element_text(size=20),
        axis.title.x= element_text(size=20),axis.title.y= element_text(size=20),
        legend.text= element_text(size=20), legend.title=element_text(size=20))

  
###DNA 测序深度 ———??? 小样，在服务器中处理
data_5_1<- read.table("D:\\backup\\our_data\\DNA\\DP\\gD1-6.DP.head_10W",header=F,sep="=")
data_5_2<- read.table("D:\\backup\\our_data\\DNA\\DP\\gD2-6.DP.head_10W",header=F,sep="=")
data_5_3<- read.table("D:\\backup\\our_data\\DNA\\DP\\gD3-6.DP.head_10W",header=F,sep="=")
data_5_4<- read.table("D:\\backup\\our_data\\DNA\\DP\\gD4-4.DP.head_10W",header=F,sep="=")
data_5_5<- read.table("D:\\backup\\our_data\\DNA\\DP\\gD5-6.DP.head_10W",header=F,sep="=")
data_5_6<- read.table("D:\\backup\\our_data\\DNA\\DP\\gD6-7-1.DP.head_10W",header=F,sep="=")
###统计DP情况
data_5_mean <- c(mean(data_5_1$V2),mean(data_5_2$V2),
                 mean(data_5_3$V2),mean(data_5_4$V2),
                 mean(data_5_5$V2),mean(data_5_5$V2))

data_5_median <- c(median(data_5_1$V2),median(data_5_2$V2),
                 median(data_5_3$V2),median(data_5_4$V2),
                 median(data_5_5$V2),median(data_5_5$V2))
#XXX<- read.table("G:\\backup\\毕业论文\\gD1-6.DP_shuf",header=F,sep="=")

print("1")

data_5_1$lable<-"gD1"
data_5_2$lable<-"gD2"
data_5_3$lable<-"gD3"
data_5_4$lable<-"gD4"
data_5_5$lable<-"gD5"
data_5_6$lable<-"gD6"
print("2")
data_5_1 <-data_5_1[,c(2,3)]
data_5_2 <-data_5_2[,c(2,3)]
data_5_3 <-data_5_3[,c(2,3)]
data_5_4 <-data_5_4[,c(2,3)]
data_5_5 <-data_5_5[,c(2,3)]
data_5_6 <-data_5_6[,c(2,3)]

gc()
print("3")
data_5<- rbind(data_5_1,data_5_2)
data_5<- rbind(data_5,data_5_3)
data_5<- rbind(data_5,data_5_4)
data_5<- rbind(data_5,data_5_5)
data_5<- rbind(data_5,data_5_6)
print("4")
rm(data_5_1)
rm(data_5_2)
rm(data_5_3)
rm(data_5_4)
rm(data_5_5)
rm(data_5_6)
gc()
print("5")
ggplot(data_5, aes(x=factor(lable), y=V2,fill=lable)) + 
  geom_violin(trim=FALSE,color="white")+stat_ydensity(bw = "nrd0")+theme_classic()+ylim(0,3)+   ###提琴???
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7"))+
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA)+   ###加入箱线图并去除离散的点
  ylab("DP")+xlab("F0 individual")+geom_hline(yintercept=10,colour="#990000")+guides(fill=FALSE)   ###添加y轴值为10 ???
data_5_p
print("6")
ggsave(data_5_p, file='DP_count.pdf', width=6, height=4)


###########   logMm 阳性基???+目标基因
library(qvalue)
library(gridExtra)
setwd("D:\\backup\\毕业论文\\logMm 计算")
list_name<-c('ENSDARG00000099517_1353864')
###  ENSDARG00000099517_1353864，，ENSDARG00000024433_61182774
i=1
file.name <- list_name[[i]]
limit <- 10000000  # 文件中存在非预期范围内的位点，剔除该位点使用的标???
p.thrd <- 0.05  # 统计检验的显著标准
## 读入 fish的genotype 预测结果
d.gene <- read.csv(file = file.name,sep = "\t",header = F)

## 根据不同的maker SNP分别输出结果???
m.pos <- unique(d.gene$V5)
print(length(m.pos))
for (m in m.pos) {
  d <- d.gene[d.gene$V5==m,]
  
  site <- unique(d$V16)
  
  ## 获取 RNA的reads 并计算log AEI
  a <- unique(d[,c(2,6,7,8,9)])
  colnames(a) <- c("id","a1","a2","an1","an2")
  a$log_AEI <- abs(log2(a$an1/a$an2))
  a$fish_id <- substr(a$id,3,4)
  
  a1 <- a[,c("fish_id","log_AEI")]

  ## 生成所有位点的phenotype ??? genotype 数据：SS为最终的理想数据格式,map为snp位置数据
  ss <- a1
  map <- data.frame(id=NA,chr=NA,pos=NA)
  map<- map[-1,]
  print(length(site))
  
  for (i in 1:length(site)) {
    s1 <- d[d$V16==site[i],c(2,4,16,19,17,22,20,24)]  ## 选择其中一个SNP位点
    s1
    colnames(s1) <- c("id","chr","pos","p1","allele1","p2","allele2","genotype")
    s1$geno2 <- ifelse(s1$allele1==s1$allele2,0,1)  ## 再次计算genotype
    s1
    
    df <- data.frame(id=paste("s",i,sep=""),chr=unique(s1$chr),pos=unique(s1$pos))
    map <- rbind(map,df)
    
    #name1 <- paste("geno",unique(s1$chr),unique(s1$pos),sep = "_")
    name1 <- paste("s",i,sep="")
    colnames(s1)[length(colnames(s1))] <-name1
    
    s1$fish_id <- substr(s1$id,3,4)
    ns1 <- ncol(s1)
    ss <- merge(ss,s1[,c(ns1,ns1-1)],by = "fish_id",all = T)
    print(i)
    
  }
  ss  # SNP 的genotype ??? phenotype文件，类似ped格式文件
  map # SNP 位置信息文件，类map格式文件
  #write.table(ss,file = paste(file.name,m,".ped.het1",sep="_"),quote = F,sep="\t",row.names = F)
  # write.table(map,file = paste(file.name,m,".map",sep="_"),quote = F,sep="\t",row.names = F)
  
  # 统计检验，获得所有位点的检验p值，并绘图（此例使用Wilcoxon rank sum test，对满足正态的数据，可是使用t test???
  ## 生成一个空dataframe
  t.d <- data.frame(chr=NA,pos=NA,t.pval=NA,w.pval=NA)
  t.d<- t.d[-1,]
  
  for (i in 1:(ncol(ss)-2)) {
    d.p <- ss[,c(1,2,i+2)]
    geo.name <-colnames(d.p)[3]
    # chr <- strsplit(geo.name,split = "_")[[1]][2] ## 获取染色体号
    # pos <- strsplit(geo.name,split = "_")[[1]][3] ## 获取SNP位置???
    chr <- map[map$id==geo.name,"chr"] ## 获取染色体号
    pos <- map[map$id==geo.name,"pos"] ## 获取SNP位置???
    
    colnames(d.p)[3] <- "genotype"
    x <- d.p$genotype
    y <- d.p$log_AEI
    
    level <- length(unique(na.omit(x)))
    
    if (level<2) {  ## 对于一??? 去除NA后只剩一种genotype的位点，需要去???
      tp <- NA
      wp <- NA
      
    }else {
      tt <- t.test(y~x,alternative="less") ## 使用单尾检验，Ha??? group 0 < group 1
      tp <- tt$p.value
      
      wt <- wilcox.test(y~x,alternative="less") # 使用单尾检验，Ha??? group 0 < group 1
      wp <- wt$p.value
    }
    
    #
    df <- data.frame(chr=chr,pos=pos,t.pval=tp,w.pval=wp)
    t.d <- rbind(t.d,df)
    
    d.p1 <- na.omit(d.p) ## 去除带有na的行
    
  }
  t.d
  t.d <- na.omit(t.d)
  #t.d <- t.d[t.d$pos>limit,]  ###  剔除位置不正确的位点
  
  # t test 结果
  t.d$t.fdr <-  p.adjust(t.d$t.pval,"BH") # global FDR矫正
  #t.d$t.qval <- qvalue(t.d$t.pval)$qvalue      # qvalue矫正
  
  # wilcox rank sum test 结果???
  t.d$w.fdr <-  p.adjust(t.d$w.pval,"BH") # global FDR矫正
  #t.d$t.qval <- qvalue(t.d$t.pval)$qvalue      # qvalue矫正
  
  min(t.d$t.fdr,na.rm = T)
  min(t.d$t.qval,na.rm = T)
  ##  t.d的是最主要的结果文???
  
  #导出文件
  
  #绘图使用wilcox fdr???(图可以使用ggsave 导出)
  t.d$col <- "fdr>0.05"
  t.d[t.d$w.fdr<=0.01,]$col <-"fdr<0.01"
  t.d[t.d$w.fdr>=0.01 & t.d$w.fdr<=0.05,]$col<-"fdr<0.05"
  mycolors<-c("#000000","#CCCCCC","#666666")
  
}
ggplot(t.d,aes(x=pos,y=-log10(t.d$w.fdr)),fill=as.factor(col))+
  geom_point(aes(color=as.factor(col),shape=as.factor(col)),size=3)+
  theme_classic()+ theme(axis.text.x = element_text(angle=45, hjust=1))+
  scale_color_manual(values = mycolors)+
  geom_hline(yintercept = -log10(p.thrd),colour="grey50")+
  ylab("-log10 FDR")+labs(title="adssl1")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="top")+###添加标题并剧???
  theme(axis.text.x = element_text(size=20),axis.text.y=element_text(size=20),
        axis.title.x= element_text(size=15),axis.title.y= element_text(size=20),legend.text= element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 10)))  ###更改图例大小

####haploview补充???
setwd("D:\\backup\\毕业论文\\haploview")
info <- read.table ("ENSDARG00000099517_1422360_1423931_.info_2",sep = "\t",header = F)
p <- read.table ("ENSDARG00000099517_1353864_.p",sep = "\t",header = T)
p_work <- p[p$pos%in%info$V2,]
colnames(info) <- c("pos_name","pos","pos_work")
info_p <- merge(p_work,info,by = "pos")
ggplot(info_p,aes(x=pos_work,y=-log10(w.pval)))+
  geom_bar(stat="identity",fill="#56B4E9",width=0.5)+
  theme_classic()+geom_hline(yintercept=-log10(0.05),colour="#CC79A7")+
  ylab("-log10(p)")+xlab("pos")

### 表达量鉴???
setwd("D:\\backup\\allele\\实验验证\\表达量比???")
data<-read.table ("adss1_pvalb4.TPM_result_work-2.txt",sep = "\t",header = T)
data_adss1 <- data[data$name=="adss1",]
data_pvalb4 <- data[data$name=="pvalb4",]

###score_sex___pvalb4
ggplot(data_pvalb4,aes(x=factor(score),y=TPM_norm.factors,fill=factor(score)))+
  geom_boxplot()+theme_classic()+
  scale_fill_manual(values=c("#CC79A7","#56B4E9", "#009E73", "#F0E442", "#D55E00"))+
  ylab("TPM")+xlab("base")+guides(fill=FALSE)+
  scale_x_discrete(labels = c("AA,TT,CC,AA,CC,AA","TA,TG,CT,AT,CT,AG","TA,TG,CT,XX,TT,XX","TT,GG,TT,XX,TT,XX","TT,GG,TT,TT,TT,GG"))   + ###改变x轴刻度名
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  theme(axis.text.x = element_text(size=20),axis.text.y=element_text(size=20),
        axis.title.x= element_text(size=20),axis.title.y= element_text(size=20))
###score  单个位点___pvalb4
ggplot(data_pvalb4,aes(x=factor(pos5),y=TPM_norm.factors,fill=factor(pos5)))+
  geom_boxplot()+theme_classic()+
  scale_fill_manual(values=c("#CC79A7","#56B4E9", "#F0E442", "#009E73","#D55E00"))+
  ylab("TPM")+xlab("base")+guides(fill=FALSE)+
  scale_x_discrete(labels = c("AA","AG","GG"))+###改变x轴刻度名
  theme(axis.text.x = element_text(size=20),axis.text.y=element_text(size=20),
        axis.title.x= element_text(size=20),axis.title.y= element_text(size=20))

###score_sex___adss1
data_adss1$Base.type<-factor(data_adss1$Base.type,levels=c("AG,AA,CC","AG,TA,TC","AA,TA,TC","GG,AA,CC","AA,AA,CC"),ordered=TRUE)
ggplot(data_adss1,aes(x=factor(Base.type),y=TPM_norm.factors,fill=factor(Base.type)))+
  geom_boxplot()+theme_classic()+
  scale_fill_manual(values=c("#CC79A7","#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7"))+
  ylab("TPM")+xlab("base")+guides(fill=FALSE)
###score  单个位点___adss1
#P1
data_adss1$pos1<-factor(data_adss1$pos1,levels=c("AG","AA","GG"),ordered=TRUE)
ggplot(data_adss1,aes(x=factor(pos1),y=TPM_norm.factors,fill=factor(pos1)))+
  geom_boxplot()+theme_classic()+
  scale_fill_manual(values=c("#CC79A7","#56B4E9", "#F0E442"))+
  ylab("TPM")+xlab("base")+guides(fill=FALSE)+
  scale_x_discrete(labels = c("AG","AA","GG"))###改变x轴刻度名
#P2
data_adss1$pos2<-factor(data_adss1$pos2,levels=c("TA","AA"),ordered=TRUE)
ggplot(data_adss1,aes(x=factor(pos2),y=TPM_norm.factors,fill=factor(pos2)))+
  geom_boxplot()+theme_classic()+
  scale_fill_manual(values=c("#CC79A7","#56B4E9"))+
  ylab("TPM")+xlab("base")+guides(fill=FALSE)+
  scale_x_discrete(labels = c("TA","AA"))###改变x轴刻度名
#P3
data_adss1$pos3<-factor(data_adss1$pos3,levels=c("TC","CC"),ordered=TRUE)
ggplot(data_adss1,aes(x=factor(pos3),y=TPM_norm.factors,fill=factor(pos3)))+
  geom_boxplot()+theme_classic()+
  scale_fill_manual(values=c("#CC79A7","#56B4E9", "#F0E442"))+
  ylab("TPM")+xlab("base")+guides(fill=FALSE)+
  scale_x_discrete(labels = c("TC","CC"))###改变x轴刻度名
###
ggplot(data_adss1,aes(x=lable_2,y=TPM_norm.factors))+
  geom_boxplot()+theme_classic()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               aes(fill=sex))+
  scale_fill_manual(values=c("#CC79A7","#56B4E9", "#F0E442", "#009E73"))+
  ylab("TPM")+xlab("adss1")


##test pvalb4
info<-read.table ("ENSDARG00000024433_61205596_61206382_.info",sep = "\t",header = F)
pvalb4.genotype<-read.table ("pvalb4.genotype",sep = "\t",header = F)
genotype.info <- merge(pvalb4.genotype,info,by="V2")
write.table(genotype.info,"pvalb4.genotype_info",row.names = F,quote = TRUE,col.names = F)



RNA_pos.all<-read.table ("D:\\backup\\our_data\\RNA\\VQSR\\10\\gene_area\\intersect_191205\\添加新的候选基因\\RNA_pos.all",sep = "\t",header = F)
XXXXXXXXXX<-unique(RNA_pos.all$V1)


####平均荧光强度
library(ggplot2)
setwd("D:\\backup\\allele\\实验验证\\注射图片")
data1 <- read.table("平均荧光强度.txt",header=T,sep = "\t")   ###单个
ggplot(data1,aes(x=factor(time),y=mean,fill=type))+
  geom_bar(position = "dodge", stat = "identity",width=0.5)+  ###分组参数
  ylab("average flourescence intensity")+xlab("time")+
  scale_fill_manual(values=c("#CC79A7","#56B4E9"))+
  theme_classic()+theme(legend.position="bottom")  ##去除背景

data2 <- read.table("多个焦面平均荧光强度.txt",header=T,sep = "\t")   ###多个
####作图前计???
mean <- aggregate(data2$light, by=list(data2$time, data2$type), FUN=mean)
sd <- aggregate(data2$light, by=list(data2$time, data2$type), FUN=sd)
len <- aggregate(data2$light, by=list(data2$time, data2$type), FUN=length)
df_res <- data.frame(mean, sd=sd$x, len=len$x)
colnames(df_res) = c("time", "type", "Mean", "Sd", "Count")
str(df_res)
df_res$Se <- df_res$Sd/sqrt(df_res$Count) ### 计算标准???

ggplot(df_res,aes(x=factor(time),y=Mean,fill=type))+
  geom_bar(position = "dodge", stat = "identity",width=0.5)+  ###分组参数
  ylab("average flourescence intensity")+xlab("time")+
  scale_fill_manual(values=c("#CC79A7","#56B4E9"))+
  geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean +Sd),position=position_dodge(.5), width=.2)+
  theme_classic()+theme(legend.position="bottom")  ##去除背景
#####  T检验_enchance
data2_24 <- data2[data2$time=="24",]
t.test(light~type,data2_24)
data2_36 <-data2[data2$time=="36",]
t.test(light~type,data2_36)
data2_48 <-data2[data2$time=="48",]
t.test(light~type,data2_48)
#####  T检验_canditate
##low-high
data<- read.table("low_high_2st_inject_2.txt",sep = "\t",header = T,row.names = 1)
colnames(data)<-c('f','nf')
fisher.test(data)
##low-high
x<-c(14,140)
y<-c(26,146)
d<-data.frame(x=x,y=y)
fisher.test(d)
##EVC-high
x<-c(14,140)
y<-c(48,5)
d<-data.frame(x=x,y=y)
fisher.test(d)
####注射后有无荧光统???
setwd("D:\\backup\\毕业论文")
data<- read.table("low_high_2st_inject.txt",sep = "\t",header = T) 
fisher.test(data)
ggplot(data,aes(x=process,y=count,fill=Fluorescence))+
  geom_bar(stat="identity",position = "fill",color="black")+ylab("freqency")+
  theme(legend.position = "top")+theme_classic()+
  theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=20),
        axis.title.x= element_text(size=15),axis.title.y= element_text(size=20))+
  scale_fill_manual(values=c("#CCCCCC","#333333"))+
  guides(shape = guide_legend(override.aes = list(size = 20)))+  ###更改图例大小
  theme(legend.position="bottom")
#####logMm模拟数据
e1<-data.frame(x=rnorm(36,mean=1.5,sd=0.5),y="heterozygous")
e2<-data.frame(x=rnorm(36,mean=0,sd=0.5),y="homozygote")
d_e<-rbind(e1,e2)
ggplot(d_e,aes(x=y,y=abs(log2(x)),fill=y))+
  geom_boxplot()+theme(legend.position="top")+theme_classic()+   ###去除背景网格+
  scale_fill_manual(values=c("#CC79A7","#56B4E9"))+
  geom_jitter(aes(fill=y),width =0.2,shape = 21,size=2)+
  ylab("Balance Score")+xlab("genotype")+guides(fill=FALSE)+
  scale_x_discrete(labels = c("het","homo"))+###改变x轴刻度名
  theme(axis.text.x = element_text(size=20),axis.text.y=element_text(size=20),
        axis.title.x= element_text(size=20),axis.title.y= element_text(size=20))


## qpcr
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






