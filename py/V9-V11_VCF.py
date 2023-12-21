import os,re
os.chdir("D:\\backup\\our_data\\RNA\\paper_data")
vcf={}
with open("SNPfisher_FLI.vcf.bed") as f1:
    for l1 in f1:
        l11=l1.rstrip("\n").split("\t")
        chr_vcf = l11[0]
        pos=l11[1]
        x1=l11[2]
        x2=l11[3]
        x3=l11[4]
        x4=l11[5]
        x5=l11[6]
        x6 = l11[7]
        x7 = l11[8]
        x8 = l11[9]
        x9 = l11[10]
        x10 = l11[11]
        y=x1+'\t'+x2+'\t'+x3+'\t'+x4+'\t'+x5+'\t'+x6+'\t'+x7+'\t'+x8+'\t'+x9+'\t'+x10
        vcf[chr_vcf,pos]=y
bed={}
with open("FLI.bed.v11") as f2:
    for l2 in f2:
        l22=l2.rstrip("\n").split("\t")
        chr_v9=l22[3]
        pos_v9=l22[4]
        chr_v11=l22[0]
        pos_v11=l22[1]
        bed[chr_v9,pos_v9]=chr_v11,pos_v11
txt=open("FLI.vcf.v11","a")
for a,b in vcf.items():
    for c,d in bed.items():
        if a==c:
            print(d,b)
            print(d, b,file=txt)
txt.close()
