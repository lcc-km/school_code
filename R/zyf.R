#2:
setwd("C:\\Users\\lcc\\Desktop")
data<-read.table("soybean.csv",sep=",",header=TRUE)
shaking<-factor(data$shaking)
light<-factor(data$light)
summary(aov(data$leaf~shaking+light))


#3:
tab <- as.table(cbind(c(271,240), c(343,216)))  #创建列联表
dimnames(tab) <- list(c("ASD", "no-ASD"),c("T", "C"))
tab
chisq.test(tab)
fisher.test(tab)


#4:
tab <- as.table(cbind(c(35,49), c(2345,28367)))  #创建列联表
dimnames(tab) <- list(c("interested", "background"),c("interested_all", "background_all"))
tab
fisher.test(tab)
