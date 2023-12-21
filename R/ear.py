import os,re
os.chdir("C:\\Users\\lcc\\Desktop\\ear")
gtf={}
#txt = open("ear_RNA-seq_result_high.txt", "a")
with open("high.txt") as f1:
    for l1 in f1:
        l11=l1.rstrip("\n").split("\t")
        E_ID = l11[0]
        name2 =l11[1]
        ncbi_ID=l11[2]
        gtf[E_ID]=name2,ncbi_ID
ear={}
with open("ear_RNA-seq_result_high.txt") as f2:
    for l2 in f2:
        l21 = l2.rstrip("\n").split("\t")
        ID=l21[0]
        name_1=l21[1]
        logFC=l21[2]
        logCPM=l21[3]
        PValue=l21[4]
        FDR=l21[5]
        ear[ID]=name_1,logFC,logCPM,PValue,FDR
for a,b in gtf.items():
    for c,d in ear.items():
        if a==c :
          print(a,b,d,file=txt)
txt.close()