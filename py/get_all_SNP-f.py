import os,re
#from collections import defaultdict
#import json
os.chdir("D:\\backup\\our_data\\DNA\\get_SNP")
#txt=("gD1-6.homozygosis","w")
#Dict={}
with open("1") as f1:
    for l1 in f1:
        l11=l1.rstrip("\n").split("\t")
        chr=l11[0]
        pos=l11[1]
        ref=l11[2]
        ale=l11[3]
        ALE=ale.split(",")
        y=len(ALE)
        DP=l11[4]
        I16=l11[5]
        i16=I16.replace("=",",").split(",")
        x1=i16[3]
        x2=i16[4]
        if ale !="<*>":
            print(chr,pos,ale,DP,x1,x2)