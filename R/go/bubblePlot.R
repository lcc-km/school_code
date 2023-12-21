######################################################
# This R script is designed for drawing Bubble plots #
# using GO/KEGG functional enrichment results from   #
# DAVID website.                                     #
#----------------------------------------------------#
# Author: < Xiao LieBao > from < Ke Yan Mao >        #
# Affiliation: Shanghai Tengyun BioTech Co.,ltd.     #
# Email: customer_service1501@tengyunbio.com         #
# Website: http://www.tengyunbio.com                 #
#                                                    #
# Date: 2019-2-16                                    #
# Version: 0.1                                       #
######################################################
#                    CAUTION                         3
#----------------------------------------------------#
# Copyright (C) 2019 by Tengyun BioTech Co.,ltd.     #
# All rights reserved.                               #
######################################################

options(stringsAsFactors = F)

# Packages. If you don't have it, install it!
library(Hmisc)
library(ggplot2)

# Functions to draw plots
DrawGOBubblePlot <- function(dat, category = "BP", top.number = 10, col="blue"){
  # Draw bubble plot using DAVID function enrichment results
  
  category = toupper(category)
  if (category == "BP"){
    main.title = "Biological Process"
  } else if (category == "CC"){
    main.title = "Cellular Components"
  } else if (category == "MF"){
    main.title = "Molecular Function"
  } else if (category == "KEGG"){
    main.title = "KEGG"
  } else {
    return("ERROR! Wrong input parameter [category].")
  }
  
  dat1 = dat[c(1:top.number),c(2,3,4,5)]
  colnames(dat1)[3] = 'GeneRatio'
  if(category == 'KEGG'){
    dat1$Term = substr(dat1$Term,10,200)
  }else{
    dat1$Term = substr(dat1$Term,12,200)
  }
  
  
  dat1$Term = capitalize(dat1$Term)
  dat1$Term = factor(dat1$Term,levels=dat1$Term[length(dat1$Term):1])
  dat1$PValue = -log10(dat1$PValue)
  
  p = ggplot(dat1,aes(GeneRatio,Term)) +
    geom_point(aes(size=Count,colour=PValue)) +
    scale_colour_gradient(low=col,high="red") + 
    labs(colour=expression(-log[10]("P Value")),size="Gene counts",  
         x="Gene Ratio",y="",title=main.title) +
    theme_bw() +
    scale_x_continuous(limits = c(0,max(dat1$GeneRatio) * 1.2)) 
  
  return(p)
}

# Read in data and generate the plots
# Biological Process
dat = read.table(file.choose(),header=T,sep="\t")
DrawGOBubblePlot(dat,"BP",10,"blue")
# Cellular Component
dat = read.table(file.choose(),header=T,sep="\t")
DrawGOBubblePlot(dat,"CC",10,"blue")
# Molecular Function
dat = read.table(file.choose(),header=T,sep="\t")
DrawGOBubblePlot(dat,"MF",10,"blue")
# KEGG Pathway
dat = read.table(file.choose(),header=T,sep="\t")
DrawGOBubblePlot(dat,"KEGG",10,"blue")

