# 17.8.11
import os
import re
os.chdir("C:/Users/49248/Desktop")
def readfasta(filename):
    fa = open('Danio_rerio.GRCz11.pep.all.fa', 'r')
    res = {}
    ID = ''
    for line in fa:
        if line.startswith('>'):
            ID = line#.strip('\n')
            res[ID] = ''
        else:
            res[ID] += line#.strip('\n')
    return res
res = readfasta('Danio_rerio.GRCz11.pep.all.fa')
regex = re.compile(r'=\d+')
uniq = {}
longest = 0
for k,v in res.items():
    i = k.split(' ',5)
    title = (i[0]).split("_")
    title = '_'.join(title[::])+'\n'
    v = [v]
    if title not in uniq:  # 注意这种生成双层字典的方法！
        uniq[title]  = v
    else:
        uniq[title] += v
max_seq = {}
# 找出一个列表中最长的字符串并只把它留下 用max
for k,v in uniq.items():
    seq = max(v, k = len)
    max_seq[k] = seq
w = open("longest.txt","w")
for k,v in max_seq.items():
    w.write(k)
    w.write(v)
w.close()