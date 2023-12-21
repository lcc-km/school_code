import os
os.chdir("D:\\backup\\our_data\\RNA\\raw_clean\\vcf\\WORK\\gene.locate\\locate")
d={}
with open ("F.union.bed.sort") as f1:
  	for l1 in f1:
          l11=l1.rstrip("\n").split("\t")
          chr =l11[0]
          locate=l11[1]
          gene_name=l11[2]
          id=l11[3]
          result=chr+":"+locate+"-"+locate
          d[chr,locate]=result
for a ,b in d.items():
        txt = open("F.union.bed.sort.work", "a")
        print(b,file=txt)
        txt.close()