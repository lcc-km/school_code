import os,re
os.chdir("D:\\backup\\our_data\\DNA\\get_SNP")
vcf={}
txt=open("gD1-6.SNP.add","w")
with open("gD1-6.vcf.work_1W") as f1:
    for l1 in f1:
        l11=l1.rstrip("\n").split("\t")
        chr=l11[0]
        pos=l11[1]
        ref=l11[2]
        ALT=l11[3]
        alt=ALT.split(",")
        y=len(alt)
        DP=l11[4]
        I16=l11[5]
        vcf[chr,pos]=ref,ALT,DP,I16,y
for a,b in vcf.items():
    if b[4]==1:
        print(a,b[0:4],file=txt)
txt.close()

