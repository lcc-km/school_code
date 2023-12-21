import os,re
import numpy as np
os.chdir("D:\\backup\\our_data\\DNA\\get_SNP\\clean.homozygosis")
gtf={}
gD={}
txt1=open("gene","a")
txt2=open("non-gene","a")
with open("Danio_rerio.GRCz11.94.gene.gtf") as f1:
    for l1 in f1:
        l11=l1.rstrip("\n").split("\t")
        chr=l11[0]
        start=l11[1]
        end=l11[2]
        pair=l11[3]
        ID=l11[4]
        name=l11[5]
        gtf[chr,start,end]=pair,ID,name
with open ("1") as f2:
    for l2 in f2:
        l22=l2.rstrip("\n").split("\t")
        CHR=l22[0]
        post=l22[1]
        snp=l22[2]
        DP=l22[3]
        read1=l22[4]
        read2=l22[5]
        gD[chr,post]=snp,DP,read1,read2
for a,b in gtf.items():
    for c,d in gD.items():
        if a[0]==c[0] and a[1] <= c[1] <= a[2]:
            print(a[0],c[1],b,d[0],d[1],file=txt1)
        else:
            print(a[0],c[1],d[0],d[1],file=txt2)
txt1.close()
txt2.close()

