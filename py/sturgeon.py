import os
os.chdir("C:\\Users\\lcc\\Desktop\\sturgeon")
from collections import OrderedDict
ID1 ={}
with open('longest_orfs.pep') as f1:
   for line in f1:
      if line.startswith('>'):
          l11 = line.rstrip("\n").split(" ")
          id =l11[0].rstrip(">")
          id = id[1::]
          length =l11[1]
          gc=l11[2]
          TRINITY=l11[3]
          ID1[id]=length,gc,TRINITY
ID2={}

with open('5.swiss.trans.high.bitscore.id') as f2:
   for l2 in f2:
       l22 = l2.rstrip("\n")
       #id2=l22[0]
       ID2[l22]=l22
txt = open("orf_high.bitscore2.txt", "a")
for a,b in ID1.items():
    for c,d in ID2.items():
        if a==c:
            print(a,b,file=txt)
txt.close()