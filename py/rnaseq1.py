import os
os.chdir("C:/Users/49248/Desktop/parallel comparisonRNA-seq/result")
from collections import OrderedDict
d1=OrderedDict()
with open('GCF_000951615.1_common_carp_genome_feature_table.txt') as f1:
    for l1 in f1:
        l11 = l1.rstrip("\n").split("\t")
        gene = l11[0]
        ID = l11[1]
        d1[gene]=ID

d2 = OrderedDict()
with open ('co_zeb_cyp_work.txt') as f2:
    for l2 in f2:
        l22=l2.rstrip("\n").split("\t")
        zebID= l22[0]
        Id =l22[1]
        if zebID not in Id:
            d2[zebID] = [Id]
        else:
             d2[zebID].append(Id)

f3 = {}
for a,b in d1.items():
    for c,d in d2.items():
        if b == c :
          f3[a] = d
v1= {}
v1 = f3

def itertransfer(v1):
  for k,value in v1.items():
      for x in value:
          yield (k,x)
for k,x in itertransfer(v1):
  print(k,x)

w = open("co_RNAseq_id_corresponding.txt", "w")
for k, v in v1.items():
   k3 =''.join(k)
   v3 ="".join(v)
   x = k3+ "\t"+v3+'\n'
   #w.write(k3)
   #w.write(v3)
   w.write(x)
w.close()