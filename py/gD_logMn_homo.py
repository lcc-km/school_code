import os,re
os.chdir("D:\\backup\\our_data\\DNA\\QC\\logMm")
SNP={}
file=['gD1-6_2','gD2-6_2','gD3-6_2','gD4-4_2','gD5-6_2','gD6-7-1_2']
for index in range(len(file)):
 with open(file[index]) as f1:
    for l1 in f1:
        l11=l1.rstrip("\n").split("\t")
        pos=l11[0]
        ref=l11[1]
        #ale=l11[2]
        #x = ale.split(",")
        #y=len(x)
        #DP=l11[3]
        z1=l11[4]
        z2=l11[5]
        #z3=l11[6]
        #z4=l11[7]
        Z=int(z1) + int(z2)
	    #index=l11[9]
        k=pos
        v=ref,Z
        SNP={}
        SNP[k]=v
        name1 = file[index] + "_homozygosis"
        txt=open(name1,"a")
        for a,b in SNP.items():
            if int(b[1]) >= 10 :
               print(a,b[0],b[1])
