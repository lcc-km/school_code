import os
os.chdir("D:\\backup\\sturgeon\\3")
blast={}
with open("5.longest_orfs.uniprot_sprot_20181228.blastp.e5.outfmt6") as f1:
    for l1 in f1:
        l11=l1.rstrip("\n").split("\t")
        ID =l11[0].split("::")
        id_orf1=ID[1]
        id_orf2 = ID[2]
        id_orf3 = ID[3]
        orf_id=id_orf1+"::"+id_orf2+"::"+id_orf3
        target=l11[1]
        bit_score=l11[11]
        V3=l11[2]
        V4=l11[3]
        V5=l11[4]
        V6=l11[5]
        V7=l11[6]
        V8=l11[7]
        V9=l11[8]
        V10=l11[9]
        V11=l11[10]
        blast[orf_id,V3,V4,V5,V6,V7,V8,V9,V10,bit_score]=target,bit_score
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
txt = open("3-1-2", "a")
for a, b in blast.items():
    for c,d in orf.items():
        if a[0]==c[2]:
            print(a[0],b)
            print(a[0],b,file=txt)
txt.close()
