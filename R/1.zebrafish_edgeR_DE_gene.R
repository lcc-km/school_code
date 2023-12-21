library(edgeR)

setwd("C:\\Users\\49248\\Desktop\\parallel comparisonRNA-seq\\28-zeb")

spe <- "zeb_28_"  ## IhVC1H=O5DNoVV

#############        readin data & keep large counts gene only      ###############
options(stringsAsFactors = FALSE)
data <-  read.csv("gene_count_28.csv", header=T, row.names=1)
keep <- rowSums(cpm(data)>2)>=2 # more than 2 counts per million in at least 2 repeats group
d <- data[keep,]
conditions <- factor(colnames(d))

### Ht4fTZVX84QyF7#,PhR*8|8DC?R;AP5DQyF7C{3F#,2N?<RTOB ##################
# ###   change colums presented names
# substring(colnames(d),first=6)
# conditions <- factor(substring(colnames(d),first=6,last=nchar(colnames(d))-1)) #get col.names and factor as conditions
# colnames(d) <- substring(colnames(d),first=6)
# conditions
#####################################################


############ <FKc1m4o2nRl gene ##########################
###  normalize by TMM and save data to exp_study
exp_study = DGEList(counts=d, group=conditions)
exp_study = calcNormFactors(exp_study,method = "TMM") #get TMM nomalization factor
exp_study$samples
exp_study

###  PCA analysis , F@9@QyF7U{Le1m4oD#J=V.<d5DO`KF3L6H
plotMDS(exp_study,method="bcv",col="firebrick")

### QyF7C;SPVX84#,HKN*IhVCBVC
BVC=0.1

###  estimating the dispersion
exp_study1 <- exp_study #get data
exp_study1= estimateCommonDisp(exp_study1)
exp_study1$common.dispersion=BVC
exp_study1 = estimateTagwiseDisp(exp_study1)
exp_study1
plotBCV(exp_study1)


######### different expression#,fisher model  ###############################

## IhVCPhR*1H=O2nRlgene5DWi1p (JJSCSZC;SPVX845DQy1>)   6~Q!R;

# opt 1: partly compare
t <- list()
for(i in 1:length(conditions)){
  t[[i]] <- c(as.character(conditions[i]),as.character(conditions[7]))
}


## opt 2: all compare
#cond <-  c("brain","gill","heart","intestine","kidney","liver","muscle")  #PhR*8|8DWi5DC{3F
#t <- list()
#time=0
#for(i in 1:(length(cond)-1)){
#  for (j in (i+1):length(cond)){
#   time=time+1
#    t[[time]] <- c(cond[i],cond[j])
 # }

#}

##set threshold  of DE gene #############################################################
thrd <- 0.05  #set threshold of difference expression
chg_rd=1 ##set threshold of log2 change fold
##set threshold#########################################################################

## 1H=O2nRlgene

for(i in 1:length(t)){
    ###  get paires #,  UbR;2=S&8y>]PhR*8|8D#!#!
    test <- t[[i]] ##get names of test pairs
    comp <- paste(test[2],test[1],sep = "-")   ## W"Rb!!1H6T=a9{N* test[2]-test[1]
    comp

    ###   exactTest & write file
    et0 <- exactTest(exp_study1,pair = test,dispersion=BVC)
    topTags(et0)
    et0.d <-data.frame(topTags(et0,n=nrow(et0$table)))
    filename1 <- paste(spe,comp,"_all_exact_FC.txt",sep ="")
    write.table(et0.d,filename1,quote=F,sep = "\t")


    # drawing DE gene plot (optioinal)
    de <- decideTestsDGE(et0,p=thrd,lfc = chg_rd)
    summary(de) ##get number of different exp gene, -1 down regulate#,1 up regulate, 0 no difference
    plotSmear(et0,de.tags = rownames(exp_study1)[as.logical(de)])
    abline(h=c(0-chg_rd,chg_rd),col="blue")

    filename <- paste(spe,comp,"_exact_FC",chg_rd,"_FDR",thrd,".pdf",sep ="")
    pdf(file = filename,width = 5,height = 5,title = comp)
    plotSmear(et0,de.tags = rownames(exp_study1)[as.logical(de)])
    abline(h=c(0-chg_rd,chg_rd),col="blue")
    title(paste(spe,comp))
    dev.off()


    #get DE gene text file
    data0 <- et0.d[abs(et0.d[,"logFC"])>chg_rd & et0.d[,"FDR"]<thrd,]
    filename <- paste(spe,comp,"_exact_FC",chg_rd,"_FDR",thrd,".txt",sep ="")
    write.table(data0,filename,quote=F,sep = "\t")

}

###1st part end here####################################
#
#
# #############    get annotation file
# anno <- read.csv("H:\\NGS_DB_software\\DB\\Zebrafish_proteomes_20150109\\uniprot-zebrafish_20150528.tab",header = T,sep="\t",fill=T)

 #get annotation and write files
# x<-match(row.names(data0),anno$Entry.name)
# x1<-anno[x,c(3,4)]
# data.a0<-data.frame(data0,x1)
