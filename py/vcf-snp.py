import os,re
os.chdir("C:/Users/lcc/Desktop/snp/zeb_normal")
with open ("SRR1524240.raw.vcf") as f:
   p = re.compile("#")
   lines = [line for line in f.readlines() if p.search(line) is None]
vcf={}
for l1 in lines:
    l11=l1.rstrip("\n").split("\t")
    chr=l11[0]
    locate=l11[1]
    REF=l11[3]
    ALT=l11[4]
    QUAL=l11[5]
    #FILTER=l11[6]
    INFO=l11[7]
    FORMAT=l11[8]
    type=l11[9]
    v=REF,ALT,QUAL,INFO,FORMAT,type
    vcf[chr,locate]=v
gene183={}
with open("183.txt") as f2:
    for l2 in f2:
        l22 = l2.rstrip("\n").split("\t")
        che183=l22[0]
        start=l22[1]
        end=l22[2]
        name=l22[3]
        gene183[che183]=start,end,name
txt=open("SRR1524240_result.txt","w")
for a,b in vcf.items():
    for c,d in gene183.items():
        if a[0]==c[0]:
            if int(d[0]) < int(a[1]) <int(d[1]):
                print(str(d[2])+"\t"+str(a[::])+"\t"+str(b[::]),file=txt)
txt.close()





