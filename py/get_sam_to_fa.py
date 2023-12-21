import os,re
os.chdir("/mnt/nas1/lucc/work/ourdata/zebarfish_gDNA/rwa_data/clean/sam/bam/full_data/SNPfisher_get_bam")
txt1 = open("gD1.fa.work", "a")
x = 1
with open("gD1.fa") as f1:
    for l1 in f1:
        if l1.startswith("@"):
            print(x, file=txt1)
            print(x)
            x += 1
        else:
              lene=l1
              lene = lene
              print(lene, file=txt1)
              #print(lene)
txt1.close()
txt2 = open("all.locate.sort.work", "w")
y = 1
with open("all.locate.sort") as f2:
    for l2 in f2:
        y = str(y)
        l22= y+"\t"+ y +"\t"+y +"\t"+ l2
        y = int(y)
        y +=1
        #print(l22)
        print(l22,file=txt2)
txt2.close()










