import os
os.chdir("C:/Users/49248/Desktop")
dict1 = []
f2=open('mart_export_work.txt', 'w')
with open ('mart_export.txt','r') as f1:
    for l1 in f1:
        l11 = l1.rstrip("\n").split('\t')
        g = l11[0]
        pv= l11[1]
        start = l11[2]
        end = l11[3]
        x = int(end)-int(start)
        z = str(x)
        y = g + '\t'+pv+'\t'+ z + '\n'
        f2.write(y)