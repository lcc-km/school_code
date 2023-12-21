import os
os.chdir("C:/Users/lcc/Desktop/meme/183/-+/low 5K+1K")
f=open('183_up.fa')
seq1={}
for line in f:
  if line.startswith('>'):
     name=line.replace('>','').split()[0]
     seq1[name]=''
  else:
     seq1[name]+=line.replace('\n','').strip()
f.close()
seq3={}
with open ("uniquen.txt") as f2:
    for l2 in f2:
       l22=l2.rstrip("\n").split("\t")
       ID=l22[0]
       chr=l22[1]
       seq3[ID]=chr
seq4={}
for a,b in seq1.items():
    for c,d in seq3.items():
        if a==d:
            seq4[c]=b
w = open("uniquen2.txt", "w")
for k, v in seq4.items():
   k += '\n'
   v += '\n'
   w.write(k)
   w.write(v)
w.close()