library(ggplot2)
library(qvalue)
library(gridExtra)
setwd("E:/OneDrive/1.research_project/#2019_zebrafish_hybrid_AEI/lucc_data_plot_try")


file.name <- "ENSDARG00000017441"
limit <- 10000000  # ?ļ??д??ڷ?Ԥ?ڷ?Χ?ڵ?λ?㣬?޳???λ??ʹ?õı?׼
p.thrd <- 0.05  # ͳ?Ƽ???????????׼

## ???? fish??genotype Ԥ??????
d.gene <- read.csv(file = file.name,sep = "\t",header = F)

## ???ݲ?ͬ??maker SNP?ֱ???????????
m.pos <- unique(d.gene$V5)

for (m in m.pos) {
  d <- d.gene[d.gene$V5==m,]
   
  site <- unique(d$V16)
  
  ## ??ȡ RNA??reads ??????log AEI
  a <- unique(d[,c(2,6,7,8,9)])
  colnames(a) <- c("id","a1","a2","an1","an2")
  a$log_AEI <- abs(log2(a$an1/a$an2))
  a$fish_id <- substr(a$id,1,2)
  
  a1 <- a[,c("fish_id","log_AEI")]
  
  # ## ѡ??genome??һ??SNPλ?? ???в???
  # s1 <- d[d$V16==site[108],c(2,4,16,19,17,22,20,24)]  ## ѡ??????һ??SNPλ??
  # s1
  # colnames(s1) <- c("id","chr","pos","p1","allele1","p2","allele2","genotype")
  # s1$geno2 <- ifelse(s1$allele1==s1$allele2,0,1)  ## ?ٴμ???genotype
  # s1
  # 
  # s1$fish_id <- substr(s1$id,1,2)
  # d.plot <- merge(a1,s1,by = "fish_id") # ??ͼ??????
  # 
  # x <- d.plot$genotype
  # y <- d.plot$log_AEI
  # 
  # test <- wilcox.test(y~x,alternative="less") # ʹ?õ?β???飬Ha?? group 0 < group 1
  # test 
  # test$p.value
  # 
  # ggplot(d.plot,aes(x=as.factor(genotype),y=log_AEI))+
  #   geom_boxplot()+
  #   geom_jitter(width = 0.2)
  # 
  # #####   ????λ?? end #########
  
  
  ## ????????λ????phenotype ?? genotype ???ݣ?SSΪ???յ????????ݸ?ʽ,mapΪsnpλ??????
  ss <- a1
  map <- data.frame(id=NA,chr=NA,pos=NA)
  map<- map[-1,]

  for (i in 1:length(site)) {
      s1 <- d[d$V16==site[i],c(2,4,16,19,17,22,20,24)]  ## ѡ??????һ??SNPλ??
      s1
      colnames(s1) <- c("id","chr","pos","p1","allele1","p2","allele2","genotype")
      s1$geno2 <- ifelse(s1$allele1==s1$allele2,0,1)  ## ?ٴμ???genotype
      s1
      
      df <- data.frame(id=paste("s",i,sep=""),chr=unique(s1$chr),pos=unique(s1$pos)) 
      map <- rbind(map,df)
      
      #name1 <- paste("geno",unique(s1$chr),unique(s1$pos),sep = "_")
      name1 <- paste("s",i,sep="")
      colnames(s1)[length(colnames(s1))] <-name1
      
      s1$fish_id <- substr(s1$id,1,2)
      ns1 <- ncol(s1)
      ss <- merge(ss,s1[,c(ns1,ns1-1)],by = "fish_id",all = T)
  
  }
  ss  # SNP ??genotype ?? phenotype?ļ???????ped??ʽ?ļ?
  map # SNP λ????Ϣ?ļ?????map??ʽ?ļ?
  write.table(ss,file = paste(file.name,m,".ped.het1",sep="_"),quote = F,sep="\t",row.names = F)
  write.table(map,file = paste(file.name,m,".map",sep="_"),quote = F,sep="\t",row.names = F)
  
  # ͳ?Ƽ??飬????????λ???ļ???pֵ??????ͼ??????ʹ??Wilcoxon rank sum test??????????̬?????ݣ?????ʹ??t test??
    ## ????һ????dataframe
    t.d <- data.frame(chr=NA,pos=NA,t.pval=NA,w.pval=NA)
    t.d<- t.d[-1,]
    
    for (i in 1:(ncol(ss)-2)) {
      d.p <- ss[,c(1,2,i+2)]
      geo.name <-colnames(d.p)[3] 
      # chr <- strsplit(geo.name,split = "_")[[1]][2] ## ??ȡȾɫ????
      # pos <- strsplit(geo.name,split = "_")[[1]][3] ## ??ȡSNPλ?ú?
      chr <- map[map$id==geo.name,"chr"] ## ??ȡȾɫ????   
      pos <- map[map$id==geo.name,"pos"] ## ??ȡSNPλ?ú?
   
      colnames(d.p)[3] <- "genotype"
      x <- d.p$genotype
      y <- d.p$log_AEI
      
      level <- length(unique(na.omit(x)))
      
      if (level<2) {  ## ????һЩ ȥ??NA??ֻʣһ??genotype??λ?㣬??Ҫȥ??
        tp <- NA
        wp <- NA
        
      }else {
        tt <- t.test(y~x,alternative="less") ## ʹ?õ?β???飬Ha?? group 0 < group 1
        tp <- tt$p.value
      
        wt <- wilcox.test(y~x,alternative="less") # ʹ?õ?β???飬Ha?? group 0 < group 1
        wp <- wt$p.value
      }
      
      #
      df <- data.frame(chr=chr,pos=pos,t.pval=tp,w.pval=wp) 
      t.d <- rbind(t.d,df)
      
      d.p1 <- na.omit(d.p) ## ȥ??????na????

    }
    t.d
    t.d <- na.omit(t.d)
    #t.d <- t.d[t.d$pos>limit,]  ###  ?޳?λ?ò???ȷ??λ??
    
    # t test ????
    t.d$t.fdr <-  p.adjust(t.d$t.pval,"BH") # global FDR????
    t.d$t.qval <- qvalue(t.d$t.pval)$qvalue      # qvalue????
    
    # wilcox rank sum test ??????
    t.d$w.fdr <-  p.adjust(t.d$w.pval,"BH") # global FDR????
    t.d$w.qval <- qvalue(t.d$w.pval)$qvalue      # qvalue????
    
    min(t.d$t.fdr,na.rm = T)
    min(t.d$t.qval,na.rm = T)
    ##  t.d????????Ҫ?Ľ????ļ?
    
    #?????ļ?
    write.table(t.d,file = paste(file.name,m,".p",sep="_"),quote = F,sep="\t",row.names = F)
  
    ## ???? log AEI??ͼ??
    plot_title <- paste(file.name,m,sep="_")
    
    p.dens <- ggplot(ss[,1:2],aes(x=log_AEI))+
      geom_density(alpha=.2, fill="#FF6666")+
      labs(title=plot_title)
    
    qqnorm(ss$log_AEI) # ??qqͼ
    p.qq <-  ggplot(ss,aes(sample=log_AEI))+
      stat_qq()+stat_qq_line()
    
    #??ͼʹ??wilcox fdr??(ͼ????ʹ??ggsave ????)
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
   
    ## ƴͼ?????ļ???
    plots <- list(p.qq,p.dens,p.pos)
    p.com <- grid.arrange(grobs=plots,layout_matrix=rbind(c(1,2),c(3,3)))
    ggsave(filename = paste(file.name,m,".pdf",sep="_"),p.com,width = 4,height = 4)
} 


