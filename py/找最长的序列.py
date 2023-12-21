import os
os.chdir("C:/Users/49248/Desktop")
from collections import OrderedDict
seq = {}
with open('Danio_rerio.GRCz11.pep.all.fa') as f1:
   for line in f1:
      if line.startswith('>'):
          name = line.replace('>', '').split()[0]
          seq[name] = ''
      else:
          seq[name] += line.replace('\n', '').strip()
          seq2={}
for k,v in seq.items():
    l = len(v)
    k += '\t'+str(l)
    seq2[k]= v
w = open("zebseq.txt", "w")
for k, v in seq2.items():
   z =k+'\t'+v+'\n'
   w.write(z)
 #  w.write(v)
w.close()

os.chdir("C:/Users/49248/Desktop")
from collections import OrderedDict

input_file = "3.txt"

dicta = OrderedDict()

with open(input_file) as ifile:
    for line in ifile:
        lista = line.rstrip("\n").split("\t")
        ensemb_id = lista[0]
        ncbi_id = lista[1]
        if ncbi_id not in dicta:
            dicta[ncbi_id] = [ensemb_id]
        else:
            dicta[ncbi_id].append(ensemb_id)
# for ke,va in dicta.items():
#	print(ke,va)
w = open("zebseqwork.txt", "w")
with open("zebseq.txt") as zeb:
    for zeb_line in zeb:
        zeb_list = zeb_line.rstrip("\n").split("\t")
        query_id = zeb_list[0]
        if query_id in dicta:
            search_ensembl_id = ";".join(dicta[query_id])
            print(search_ensembl_id + "\t" + zeb_line.rstrip("\n"))
            w.write(search_ensembl_id + "\t" + zeb_line.rstrip("\n")+'\n')
w.close()