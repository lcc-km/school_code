BiocManager::install("pheatmap")
library(pheatmap) 
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")
head(test)
pheatmap(test)
pheatmap(test, scale = "row",display_numbers=T,fontsize=15
         ,treeheight_row=100,treeheight_col=20)


result <- pheatmap(test)
summary(result)
