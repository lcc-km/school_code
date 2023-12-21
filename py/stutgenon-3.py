import os
os.chdir("C:\\Users\\lcc\\Desktop\\sturgeon\\2")
orf_high={}
with open ("orf_high.bitscore_un_com.txt") as f1:
    for l1 in f1:
        l11=l1.rstrip("\n").split("\t")
        x=l11[0]
        y=x.split("::")
        id=y[1]
        id_g=y[2]
        id_m=y[3]
        type=l11[1]
        len=l11[2]
        orf_high[id,id_g,id_m]=type,len
long_orf={}
with open ('longest_orfs.pep') as f2:
    for l2 in f2:
        if l2.startswith('>'):
            l22 = l2.rstrip("\n").split("::")
            id2=l22[1]
            id_g2=l22[2]
            x2=l22[3]
            y2=x2.split(" ")
            id_m2=y2[0]
            len2=y2[2]
            tpye2=y2[1]
            long_orf[id2,id_g2,id_m2]=tpye2,len2
txt = open("orf_high_uncom&longest2.txt", "a")
for a,b in orf_high.items():
    for c,d in long_orf.items():
        if a[0]==c[0] and d[0]== 'type:complete':
            print(a,b,c,d)
            print(a, b, c, d, file=txt)
txt.close()



