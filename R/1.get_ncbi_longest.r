### change input folder & filename here
setwd("D:\\Qsync_folder\\manuscript_conference\\2017_carp\\get_longest_protein")
fea <- read.csv("1.ncbi_table_input.txt",head=T,sep="\t")

### get protein only 
protein <- fea[fea$class=="with_protein",]

## get unique gene names
ugid <- unique(protein$GeneID)

### get longest protein name, generate a matrix
keep <- NULL
for (i in 1:length(ugid)) {
  len <- max(protein[protein$GeneID==ugid[i],]$product_length)
  pep <- protein[protein$GeneID==ugid[i] & protein$product_length==len,]
  geneid <- paste("gid_",as.character(pep[1,]$GeneID),sep = "")
  product <- as.character(pep[1,]$product_accession)
  keep <- rbind(keep,matrix(c(product,geneid),nrow = 1))
}

write.table(keep,file="2.longest_pep_id.txt",sep = "\t",quote = F,row.names = F,col.names = F)




