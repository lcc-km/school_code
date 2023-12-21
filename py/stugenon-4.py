import os
os.chdir("D:\\backup\\sturgeon\\3")
trans_high_bitscore={}
with open ("5.swiss.trans.high.bitscore") as f1:
    for l1 in f1:
        l11=l1.rstrip("\n").split("\t")
        trans_id=l11[1]
        gene=l11[2]
        pro_name=l11[3]
        gene_name=l11[4]
        V12=l11[5]
        V2=l11[6]
        merge=l11[7]
        trans_high_bitscore[trans_id]=V12,V2
orf={}
with open ('orf_high_uncom&longest2_work.txt') as f2:
    for l2 in f2:
        l22=l2.rstrip("\n").split("\t")
        id2=l22[0]
        type=l22[1]
        id3=l22[2]
        idx=id3.split("::")
        id33=idx[0]
        type2=l22[3]
        orf[id2,type,id3,type2,id33]=id3,id33
txt = open("3-2-2", "a")
for a,b in trans_high_bitscore.items():
    for c,d in orf.items():
        if a==d[1]:
            print(d[1],b)
            print(d[1],b, file=txt)
txt.close()


