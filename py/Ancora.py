import os,re
os.chdir("D:/backup/Ancora")
all={}
with open ("HCNE_danRer11_all.bed") as f1:
    for l1 in f1:
        l11=l1.rstrip("\n").split("\t")
        chr_all=l11[0]
        start_all=l11[1]
        end_all=l11[2]
        all[chr_all,start_all]=start_all,end_all
d183={}
with open("183_5000.bed") as f2:
    for l2 in f2:
        l22=l2.rstrip("\n").split("\t")
        chr_183=l22[0]
        start_183=l22[1]
        end_183=l22[2]
        gene_id=l22[3]
        d183[chr_183,start_183]=start_183,end_183,gene_id
txt=open("all$183_5000","w")
for a,b in all.items():
    for c,d in d183.items():
        if a[0]==c[0]:
            #if int(d[1]) - int(d[0]) -20000 < 10000:
              if int(d[0]) < int(b[0]) < int(b[1]) < int(d[1]) :
                print(str(a) + str(b[1])+d[2],file=txt)
txt.close()