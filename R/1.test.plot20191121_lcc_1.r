library(ggplot2)
library(qvalue)
library(gridExtra)
setwd("D:\\backup\\our_data\\RNA\\VQSR\\10\\gene_area\\intersect_191205\\result\\try")

list_name<-c("ENSDARG00000009342")
###  ENSDARG00000068214 ,'ENSDARG00000081758','ENSDARG00000095507','ENSDARG00000100753','ENSDARG00000116806'
i=1
for (i in 1:length(list_name)) 
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

  # ## 选择genome上一个SNP位点 进行测试
  # s1 <- d[d$V16==site[108],c(2,4,16,19,17,22,20,24)]  ## 选择其中一个SNP位点
  # s1
  # colnames(s1) <- c("id","chr","pos","p1","allele1","p2","allele2","genotype")
  # s1$geno2 <- ifelse(s1$allele1==s1$allele2,0,1)  ## 再次计算genotype
  # s1
  #
  # s1$fish_id <- substr(s1$id,1,2)
  # d.plot <- merge(a1,s1,by = "fish_id") # 绘图用数???
  #
  # x <- d.plot$genotype
  # y <- d.plot$log_AEI
  #
  # test <- wilcox.test(y~x,alternative="less") # 使用单尾检验，Ha??? group 0 < group 1
  # test
  # test$p.value
  #
  # ggplot(d.plot,aes(x=as.factor(genotype),y=log_AEI))+
  #   geom_boxplot()+
  #   geom_jitter(width = 0.2)
  #
  # #####   单个位点 end #########


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
  write.table(ss,file = paste(file.name,m,".ped.het1",sep="_"),quote = F,sep="\t",row.names = F)
  write.table(map,file = paste(file.name,m,".map",sep="_"),quote = F,sep="\t",row.names = F)

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
    write.table(t.d,file = paste(file.name,m,".p",sep="_"),quote = F,sep="\t",row.names = F)

    ## 针对 log AEI绘图???
    plot_title <- paste(file.name,m,sep="_")

    p.dens <- ggplot(ss[,1:2],aes(x=log_AEI))+
      geom_density(alpha=.2, fill="#FF6666")+
      labs(title=plot_title)

    qqnorm(ss$log_AEI) # 绘qq???
    p.qq <-  ggplot(ss,aes(sample=log_AEI))+
      stat_qq()+stat_qq_line()

    #绘图使用wilcox fdr???(图可以使用ggsave 导出)
    t.d$col <- "black"
    if (min(t.d$w.fdr)<=p.thrd) {
       t.d[t.d$w.fdr<=p.thrd,]$col <-"red"
    }

    mycolors<-c("black","red")

    p.pos <-
      ggplot(t.d,aes(x=pos,y=-log10(t.d$w.fdr)))+
      geom_point(aes(color=as.factor(col)),shape=1)+
      theme_classic()+
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
      scale_color_manual(values = mycolors)+
      geom_hline(yintercept = -log10(p.thrd),colour="grey50")+
      ylab("-log10 FDR")+labs(title=plot_title)

    ## 拼图导出文件???
    plots <- list(p.qq,p.dens,p.pos)
    p.com <- grid.arrange(grobs=plots,layout_matrix=rbind(c(1,2),c(3,3)))
    ggsave(filename = paste(file.name,m,".pdf",sep="_"),p.com,width = 4,height = 4)
}
print("finish")


print("finish all")

