setwd("C:\\Users\\lcc\\Desktop")
library(ggpubr)
library(ggplot2)
data <- read.csv("Anisomycin 19-4-10_A.txt",header=TRUE,sep = "\t")

background <- data[which(data$Sample.Name == "ELF2"),]    #### background?????????
background_mock <- background[which(background$Target.Name == "AW"),] #### mock???
background_handle <- background[which(background$Target.Name == "AM"),] #### handle???
Ct_mock_mean <- mean(background_mock$C??)                 ### mock???ct????????????
Ct_handle_mean <- mean(background_handle$C??)            ### handle???ct????????????

datawork <- data[,c(2,3,10)]
datawork_mock <- datawork[which(datawork$Target.Name =="AW"),]  ### ??????????????????mock??????
datawork_handle <- datawork[which(datawork$Target.Name =="AM"),]   ###  ??????????????????handle??????
data_mock_delta <- data.frame(datawork_mock$Sample.Name,datawork_mock$C??-Ct_mock_mean)
data_handle_delta <- data.frame(datawork_handle$Sample.Name,datawork_handle$C??-Ct_handle_mean)
names(data_mock_delta)<-c("Sample.Name","C??")
names(data_handle_delta)<-c("Sample.Name","C??")

##### deltaC?? C??  #####
##for(i in 1:nrow(data_handle_delta)){ t.test(data_handle_delta$C??,mu=Ct_handle_mean)


#compare_means( C??~ Sample.Name , data = data_handle_delta , method = "t.test",
#               ref.group= "ELF2" )
#p <- ggboxplot(data_handle_delta, x="Sample.Name", y = "C??",
#               palette = "npg")
#p+stat_compare_means(method = "t.test",aes(label = "p.signif"))



########
###C??_mean <-data.frame(tapply(data_mock_delta$C??,data_mock_delta$Sample.Name,mean,simplify = FALSE) )
###names(C??_mean) <-("mean")


#####   delta_delta_C??   #####
(
  if(data_mock_delta$Sample.Name == data_handle_delta$Sample.Name)
    data_delta_delta_C?? <- data.frame(data_handle_delta$Sample.Name,
                                      -(data_handle_delta$C?? - data_mock_delta$C??) )
)

names(data_delta_delta_C??)<-c("Sample.Name","C??")
compare_means( C??~ Sample.Name , data = data_delta_delta_C?? , method = "t.test",
               ref.group= "ELF2" )

p <- ggboxplot(data_delta_delta_C??, x="Sample.Name", y = "C??",
               palette = "npg", width = 0.5,size = 0.8,
               font.y=20,font.xtickslab = 20,font.ytickslab = 20,
               color ="Sample.Name",
               add = "jitter",
               ylab = "????C??",xlab ="Gene name")
for (i in 1:length(data_delta_delta_C??$Sample.Name)){

p+stat_compare_means(method = "t.test",aes(label = ..p.format..),
                     comparisons = c(i=4,i)
)
}
p+theme(legend.position='none') + labs(title = "Comparison of different expression level")
