import os
os.chdir("C:/Users/49248/Desktop")
from collections import OrderedDict
with open('longest') as f1:
    seq = {}
    for line in f1:
        if line.startswith('>'):
            name = line.split()[0]
            seq[name] = ''
        else:
            seq[name] += line.strip()

input_file="mart_export.txt"
dicta=OrderedDict()
with open(input_file) as ifile:
    for line in ifile:
          lista=line.rstrip("\n").split("\t")
          ensemb_id=lista[0]
          ncbi_id=lista[1]
          dicta[ncbi_id] = ensemb_id

seq2 = {}
for a,b in list(seq.items()):  #加上了一个list，是因为之前报错dictionary changed size during iteration
    for c,d in dicta.items():
        if a == c:
            seq2[d]  = b

write_file = os.getcwd() + '\dit.txt'
print (write_file)
output = open(write_file, 'w')

for i in seq2:
    print(i, seq2[i])
    write_str = str(i) + '\n'+ seq2[i]+'\n'
    output.write(write_str)
output.close()

