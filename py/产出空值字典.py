import os,re
os.chdir("C:/Users/lcc/Desktop/meme/183/-+/5k+1k/work/end")
seq={}
with open("background183_up_random.fa") as f1:
   for line in f1:
      if line.startswith('>'):
          name = line.split()[0]
          seq[name] = ''
      else:
          seq[name] += line.replace('\n', '').strip()
          seq2={}
for key in list(seq.keys()):
   if not seq.get(key):
      del seq[key]
w = open("background183_up_random_work.fa", "w")
for k, v in seq.items():
   k+='\n'
   v += '\n'
   w.write(k)
   w.write(v)
w.close()




















