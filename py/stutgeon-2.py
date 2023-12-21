import os
os.chdir("C:\\Users\\lcc\\Desktop\\sturgeon")
from collections import OrderedDict
ID={}
with open('longest_orfs.pep') as f1:
   for l1 in f1:
       if l1.startswith('>'):
           l11 = l1.rstrip("\n").split("::")
           l12 = l11[3].split(" ")
           id = l11[1]
           id_m = l12[0]
           len = l12[2]
           type = l12[1]
           ID[id,id_m]=type,len
txt = open("same_id_type_complete&other.txt", "a")
for a, b in ID.items():
    for c,d in ID.items():
     if a[0]==a[0] and b[0] != 'type:complete':
         if a[0] ==c[0] and d[0] == 'type:complete':
            print(a,b,c,d,file=txt)
txt.close()