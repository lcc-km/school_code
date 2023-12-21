import os
os.chdir("C:/Users/49248/Desktop")
f = open('1.txt','r')
line = ''
line = f.readlines()
for i in line:
    i = ''.join(line)
    i = i.split()
ID1 = i[0::12]
ID2 = i[1::12]
e = i[10::12]
result = []
for (a,b,c) in zip(ID1,ID2,e):
    result.append(a)
    result.append(b)
    result.append(c)
    result += '\n'
result2 = "      ".join(result)
w = open("ruesult.txt","w")
w.write(result2)
w.close()