import os,re
import gc
os.chdir("D:\\backup\\allele\\chip-seq\\download\\bed")
vcf={}
with open("GSM915196_H3K4me1_48hpf_Zv9_reads.bed") as f1:
    for l1 in f1:
        l11 = l1.rstrip("\n").split("\t")
        chr = l11[0]
        start = l11[1]
        end= l11[2]
        x1=l11[3]
        x2=l11[4]
        x3=l11[5]
        vcf[chr,start,end,x1,x2,x3]=chr,start,end