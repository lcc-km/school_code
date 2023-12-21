import os,re
from itertools import groupby
os.chdir("D:\\backup\\our_data\\DNA\\get_SNP")
from collections import OrderedDict
seq = {}
D={}
with open('gD1.work.work.work_2') as f1:
   for line in f1:
      if line.startswith('>'):
          name = line.replace('>','').split()[0]
          seq[name] = ''
      else:
          seq[name] += line.strip()
for k,v in seq.items():
    kk=k.split("_")
    seat = kk[2]
    seat = int(seat)
    vv=v.split("%")
    length = len(vv)
    length = length-1
    for i in range(length):
        l12=vv[i]
        l12=l12.split("\t")
        chr = l12[0]
        locus = int(l12[1])
        quality = l12[2]
        read = l12[3]
        ss = [''.join(list(g)) for k, g in groupby(quality, key=lambda x: x.isdigit())]
        ss_0 = int(ss[0])
        if ss[1] == "M" and locus + ss_0 >= seat:
             choose = seat - locus
             SNP = read[choose]
             D[chr, seat,locus,i] = SNP