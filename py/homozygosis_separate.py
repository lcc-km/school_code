import os,re
os.chdir("D:\\backup\\our_data\\DNA\\QC")
SNP={}
with open("gD1-6.SNP_intermediate_file") as f1:
    for l1 in f1:
        l11=l1.rstrip("\n").split("\t")
        chr=l11[0]
        pos=l11[1]
        ref=l11[2]
        ale=l11[3]
        x = ale.split(",")
        y=len(x)
        DP=l11[4]
        z1=l11[5]
        z2=l11[6]
        z3=l11[7]
        z4=l11[8]
        k=chr,pos
        v=ref,ale,DP,z1,z2,z3,z4,x,y
        SNP[k]=v
print("read finish")
txt1=open("gD1-6.SNP_homozygosis","w")
txt2=open("gD1-6.SNP_discard","w")
for a,b in SNP.items():
    if b[8] == 2:
        print(a,b[0:7],file=txt1)
    else:
        print(a,b[0:7],file=txt2)
txt1.close()
txt2.close()
